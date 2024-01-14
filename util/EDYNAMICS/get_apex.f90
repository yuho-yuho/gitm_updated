!-----------------------------------------------------------------------
      subroutine apxparm(date_run)
! 
     use fieldline_s_module,only: fieldline_s1,fline_s1, &
         fieldline_s2,fline_s2
     use fieldline_r_module,only: fieldline_r,fline_r
     use fieldline_p_module,only: fieldline_p,fline_p
     use qd_module,only: qd_grid,qd,lat_qd_mp,lon_qd_mp,lon_qd_ed, &
        hgt_qd_mp,hgt_qd_ed,lat_qd_ed,lat_qd_mp,g13,g23
     use params_module, only: ylatm,hgt_fix,nhgt_fix, &
        ylonm,nmlon,nmlat_h,nmlatS2_h,rtd,h0,r0,re,rho,rho_s, &
	nlat_qd,hgt_fix_r,pi,nhgt_fix_r,nglon,nglat,glon,glat
     use apex,only: apex_mka, apex_mall, apex_q2g,ggrid
     use delB_module,only: bhat_g,qdlon_g,qdlat_g,h_leo
 	
     implicit none	
! 
     real,intent(in) :: date_run 
!        
! to set up interpolation arrays     
!     integer, parameter:: nalt=5,nlatp1= 90, nlonp2= 180   
!     integer, parameter:: nvert
!     real, parameter :: dellat = 180./float(nlatp1-1),& 
!                        dellon = 360./float(nlonp2-2)
!     real :: gplat(nlatp1),gplon(nlonp2),gpalt(nalt)
!  
! for generating interpolation grid in apex 
     integer :: nlat,nlon,nalt  
!      integer, parameter:: nvert = 50, &
!        mxlat= 180,mxlon=360,mxalt=200
     integer, parameter:: nvert = 40, &
       mxlat= 121,mxlon=201,mxalt=8
     real, parameter :: &
        glatmin=-90.,glatmax=90.,  &
        glonmin=-180.,glonmax=180.,&
	altmin=0.,altmax=1200.
     real :: gplat(mxlat),gplon(mxlon),gpalt(mxalt)  
!
     integer :: i,j,k,ist,it,jt,itp,jtp,itm,jtm,i_start,j_start,&
        nmlath_in,ilat,isn 
     real ::  qdlat,qdlon,alt,gdlat,gdlon,hr,ds_on_fldline,fldlon,&
        fldlat,dlat,dlon,wgtp,wgtm,rho_star,sqrt_frho_star,fac1,fac1_star, &
	fac_r,dlatm
     logical found_lon,found_lat
     
     logical, parameter :: debug =.true.

! scalar arguments returned by APXMALL:    
     real :: bmag,si,alon,xlatm,vmp,w,d,be3,sim,xlatqd,f, &
       sqrt_frho,dsi,wm
! non-scalar arguments returned by APXMALL:
     real :: &
       b(3),bhat(3), &
       d1(3),d2(3),d3(3), & 
       e1(3),e2(3),e3(3),  &
       g1(3),g2(3),g3(3),  &
       f1(3),f2(3),f3(3),gradphi(3),gradrho(3),kxvec(3), &
       am1(3),am2(3),am3(3),dlonm
     real, parameter :: kunit(3) = (/0,0,1/) ! check with Art direction of unit vector
     
     integer :: itmp,jj
    
!
      dlonm = 2.*pi/float(nmlon) 
!
! Specify grid values for interpolating arrays SUBROUTINE APXMKA
!
!      gpalt=(/85,110,140,170,200,1000/)
!      gpalt=(/75,110,170,400,5000/)
!      do j=1,nlatp1
!     	gplat(j) = (j-1)*dellat - 90.
!      enddo
!      do i=1,nlonp2
!     	gplon(i) = (float(i)-1.5)*dellon - 180.
!      enddo

!   Make a global lat,lon,alt grid for use in later calls (optional)
      call ggrid(nvert,glatmin,glatmax,glonmin,glonmax,altmin,altmax, &
                 gplat,gplon,gpalt,mxlat,mxlon,mxalt,nlat,nlon,nalt)
      if(debug) write(6,*) 'done ggrid'
!
!  Initialize interpolation arrays, but do not write them
!     SUBROUTINE APXMKA (MSGUN, EPOCH, GPLAT,GPLON,GPALT,NLAT,NLON,NALT,
!    +                  WK,LWK, IST)
!            DIMENSION GPLAT(NLAT), GPLON(NLON), GPALT(NALT), WK(LWK)
!            DIMENSION GPLAT(*),GPLON(*),GPALT(*), EPOCH(*), WK(*),
!
      call apex_mka (date_run, gplat,gplon,gpalt,nlat,nlon,nalt,ist)
      if (ist /= 0) call shutdown('apxmka')
      
! 
!  Convert from quasi-dipole to geodetic coordinates, APXQ2G
!          (input magnetic, output geodetic) is the functional inverse
!          of APXALL or APXMALL (input geodetic, output magnetic). 
     do isn = 1,2
      do i=1,nmlon
        do j=1,nmlatS2_h 
	   do k=1,fline_s2(i,j,isn)%npts   
	     qdlat = fline_s2(i,j,isn)%mlat_qd(k)*rtd    ! get quasi dipole latitude 
	     qdlon = fline_s2(i,j,isn)%mlon_qd(k)*rtd    ! get quasi dipole longitude 
	     alt =fline_s2(i,j,isn)%hgt_pt(k)*1.e-3      ! height convert from [m] to [km]
             call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
             if (ist /= 0)  then
	        write(6,*) 's2 qdlon/lat/alt= ',qdlon,qdlat,alt
	        stop 'apxq2g'
	     end if	
!	     
	     fline_s2(i,j,isn)%glon(k) = gdlon  ! save geog. longitude [deg]
	     fline_s2(i,j,isn)%glat(k) = gdlat  ! save geog. latitude  [deg]
           enddo
         enddo ! end latitude loop
	 
	   
        do j=1,nmlat_h 
	   do k=1,fline_s1(i,j,isn)%npts   
	     qdlat = fline_s1(i,j,isn)%mlat_qd(k)*rtd    ! get quasi dipole latitude 
	     qdlon = fline_s1(i,j,isn)%mlon_qd(k)*rtd    ! get quasi dipole longitude 
	     alt =fline_s1(i,j,isn)%hgt_pt(k)*1.e-3      ! height convert from [m] to [km]
             call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
             if (ist /= 0)  then
	        write(6,*) 's1 qdlon/lat/alt= ',qdlon,qdlat,alt
	        stop 'apxq2g'
	     end if	
	    
	     fline_s1(i,j,isn)%glon(k) = gdlon  ! save geog. longitude [deg]
	     fline_s1(i,j,isn)%glat(k) = gdlat  ! save geog. latitude  [deg]
           enddo
	   
	   do k=1,fline_r(i,j,isn)%npts   
	     qdlat = fline_r(i,j,isn)%mlat_qd(k)*rtd    ! get quasi dipole latitude 
	     qdlon = fline_r(i,j,isn)%mlon_qd(k)*rtd    ! get quasi dipole longitude 
	     alt =fline_r(i,j,isn)%hgt_pt(k)*1.e-3      ! height convert from [m] to [km]
             call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
             if (ist /= 0)  then
	        write(6,*) 'r qdlon/lat/alt= ',qdlon,qdlat,alt
	        stop 'apxq2g'
	     end if	
!	     
	     fline_r(i,j,isn)%glon(k) = gdlon  ! save geog. longitude [deg]
	     fline_r(i,j,isn)%glat(k) = gdlat  ! save geog. latitude  [deg]
           enddo
	   
	   do k=1,fline_p(i,j,isn)%npts   
	     qdlat = fline_p(i,j,isn)%mlat_qd(k)*rtd    ! get quasi dipole latitude 
	     qdlon = fline_p(i,j,isn)%mlon_qd(k)*rtd    ! get quasi dipole longitude 
	     alt =fline_p(i,j,isn)%hgt_pt(k)*1.e-3      ! height convert from [m] to [km]
             call apex_q2g(qdlat,qdlon,alt,gdlat,gdlon,ist)
             if (ist /= 0)  then
	        write(6,*) 'r qdlon/lat/alt= ',qdlon,qdlat,alt
	        stop 'apxq2g'
	     end if	
!	     
	     fline_p(i,j,isn)%glon(k) = gdlon  ! save geog. longitude [deg]
	     fline_p(i,j,isn)%glat(k) = gdlat  ! save geog. latitude  [deg]
           enddo
         enddo ! end latitude loop
	 
       enddo
      enddo
      if(debug) write(6,*) 'done apxq2g'
! 
! now the coordinate system is [phi_m,rho=coslam_m,+/-h]
! the base vectors are d' and e'
! so there are new base vectors di' and ei'
! with d1' = d1; d2' = d2; e3' = e3; D; B0; 
! what is different is d3'; e1'; e2'
! d3' = -k^/(DsinI)
! e1' = d2'xd3' = (R k^ x grad(rho))/(sqrt(1-3/4rho^2)D sinI)
! e2' = d3'xd1' = (R rho grad(rho) x k^)/(D sinI)  
! Volume and area vectors:
! Wm'  = |k^ dot grad(phi_m) x grad(rho)|^-1 = R^2 rho/(D |sinI| sqrt(1-3/4rho^2))
! am1' = Wm' grad(phi_m)= R d1' /(D |sinI| sqrt(1-3/4rho^2))
! am2' = Wm' grad(rho_m)= R rho d2' /(D |sinI|)
! am3' = -/+Wm' k^      = R^2 rho d3' /(sqrt(1-3/4rho^2)) 
! with R = r0 =Re + h0 =Re + hr
! 
     do isn = 1,2
       do i = 1,nmlon
!      
         do j = 1,nmlatS2_h
!	 
	   do k=1,fline_s2(i,j,isn)%npts 
	     gdlon = fline_s2(i,j,isn)%glon(k)  ! geog.longitude should be in [deg]
	     gdlat = fline_s2(i,j,isn)%glat(k)  ! geog. latitude should be in [deg]
	  
             hr = h0*1.e-3                      ! reference height convert from [m] to [km]
	     alt =fline_s2(i,j,isn)%hgt_pt(k)*1.e-3  ! height convert from [m] to [km]
!	     if(alt < hr ) then
!	        write(6,*) 'hgt ',hr,alt,i,j,k
!	     endif	
             call apex_mall ( &
     	        gdlat,gdlon,alt,hr,  &		   !Inputs
     	        b,bhat,bmag,si,         &		   !Mag Fld
     	        alon,		      & 		   !Apx Lon
     	        xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	        xlatqd,f,f1,f2 ,f3,g1,g2,g3, ist)			   !Qsi-Dpl
              if (ist .ne. 0) stop 'apex ist'
!	      
! these are the same using the "new" coordinate system and the one from the paper
	      fline_s2(i,j,isn)%Vmp(k)  = vmp         ! magnitude potential Tm (diagnostic for ds calculation)
	      fline_s2(i,j,isn)%Bmag(k) = bmag*1.e-9  ! magnitude of magnetic field, convert from nT to T
	      fline_s2(i,j,isn)%sinI(k) = si          ! sin(I)
	      fline_s2(i,j,isn)%bo(:,k) = b*1.e-9     ! magnetic field components (east, north, up), up positive [T]
	      fline_s2(i,j,isn)%be3(k)  = be3*1.e-9   ! B0= Be3*e3, convert from nT to T 
	      fline_s2(i,j,isn)%D(k)    = d           ! D 
	      fline_s2(i,j,isn)%F(k)	= f	      ! F
	      fline_s2(i,j,isn)%d1(:,k) = d1          ! components (east, north, up) of base vectors
	      fline_s2(i,j,isn)%d2(:,k) = d2          ! components (east, north, up) of base vectors
	      fline_s2(i,j,isn)%e3(:,k) = e3          ! components (east, north, up) of base vectors
	      fline_s2(i,j,isn)%d1d2(k) = dot_product(d1,d2)
	      fline_s2(i,j,isn)%d2d2(k) = dot_product(d2,d2)
	      fline_s2(i,j,isn)%d1d1(k) = dot_product(d1,d1)  ! diagnostic
	      fline_s2(i,j,isn)%e1g1(k) = dot_product(e1,g1)  ! diagnostic
	      fline_s2(i,j,isn)%e2g1(k) = dot_product(e2,g1)  ! diagnostic
	      fline_s2(i,j,isn)%e2g2(k) = dot_product(e2,g2) 
	      b=  b*1.e-9
	      fline_s2(i,j,isn)%bg1(k) = dot_product(b,g1)/fline_s2(i,j,isn)%Bmag(k)  ! diagnostic
	      fline_s2(i,j,isn)%bg2(k) = dot_product(b,g2)/fline_s2(i,j,isn)%Bmag(k)  ! diagnostic
	      fline_s2(i,j,isn)%e2k(k)  = e2(3)  ! k unit upward vector  
	      fline_s2(i,j,isn)%e1g2(k) = dot_product(e1,g2)
	      fline_s2(i,j,isn)%e1k(k)  = e1(3)  ! k unit upward vector
!	      write(22,'(2(x,i4),5(x,f12.5))') i,j,gdlon,gdlat,alt,fline_s2(i,j,isn)%d1d2(k),d
!	      write(22,'(2(x,i4),5(x,f12.5))') i,j,gdlon,gdlat,alt,fline_s2(i,j,isn)%mlon_qd(k),&
!	       fline_s2(i,j,isn)%mlat_qd(k)
!	       
	    enddo ! end of height loop
	  enddo ! end of latitude loop
!      
         do j = 1,nmlat_h
!	 
	   do k=1,fline_s1(i,j,isn)%npts 
	     gdlon = fline_s1(i,j,isn)%glon(k)  ! geog.longitude should be in [deg]
	     gdlat = fline_s1(i,j,isn)%glat(k)  ! geog. latitude should be in [deg]
	   
             hr = h0*1.e-3		      ! reference height convert from [m] to [km]
	     alt =fline_s1(i,j,isn)%hgt_pt(k)*1.e-3  ! height convert from [m] to [km]
!	     if(alt < hr ) then
!	   	write(6,*) 'hgt ',hr,alt,i,j,k
!	     endif	
             call apex_mall ( &
     	   	gdlat,gdlon,alt,hr,   &		          !Inputs
     	   	b,bhat,bmag,si, 	&		  !Mag Fld
     	   	alon,		      & 		  !Apx Lon
     	   	xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	   	xlatqd,f,f1,f2,f3,g1,g2,g3 , ist)			 !Qsi-Dpl
              if (ist .ne. 0) stop 'apex ist'
	      
	      fline_s1(i,j,isn)%Vmp(k)  = vmp         ! magnitude potential Tm (diagnostic for ds calculation)
	      fline_s1(i,j,isn)%Bmag(k) = bmag*1.e-9  ! magnitude of magnetic field, convert from nT to T
	      fline_s1(i,j,isn)%sinI(k) = si          ! sin(I)
	      fline_s1(i,j,isn)%bo(:,k)   = b*1.e-9     ! magnetic field components (east, north, up), up positive [T]
	      fline_s1(i,j,isn)%be3(k)  = be3*1.e-9   ! B0= Be3*e3, convert from nT to T 
	      fline_s1(i,j,isn)%D(k)	= d	      ! D 
	      fline_s1(i,j,isn)%F(k)	= f	      ! F
	      fline_s1(i,j,isn)%d1(:,k) = d1	      ! components (east, north, up) of base vectors
	      fline_s1(i,j,isn)%d2(:,k) = d2	      ! components (east, north, up) of base vectors
	      fline_s1(i,j,isn)%e3(:,k) = e3	      ! components (east, north, up) of base vectors
	      fline_s1(i,j,isn)%d1d2(k) = dot_product(d1,d2)
	      fline_s1(i,j,isn)%d1d1(k) = dot_product(d1,d1)
	      fline_s1(i,j,isn)%d2d2(k) = dot_product(d2,d2)
	      fline_s1(i,j,isn)%e2g2(k) = dot_product(e2,g2)
	      fline_s1(i,j,isn)%e2k(k)  = e2(3)  ! k unit uoward vector
	      fline_s1(i,j,isn)%e1g2(k) = dot_product(e1,g2)
	      fline_s1(i,j,isn)%e1k(k)  = e1(3)  ! k unit uoward vector
	      
           enddo ! end of height loop
!
	   do k=1,fline_r(i,j,isn)%npts ! also for lowest level since needed for Jr calculation
	     gdlon = fline_r(i,j,isn)%glon(k)  ! geog.longitude should be in [deg]
	     gdlat = fline_r(i,j,isn)%glat(k)  ! geog. latitude should be in [deg]
	   
             hr = h0*1.e-3		      ! reference height convert from [m] to [km]
	     alt =fline_r(i,j,isn)%hgt_pt(k)*1.e-3  ! height convert from [m] to [km]
!	     if(alt < hr ) then
!	   	write(6,*) 'hgt ',hr,alt,i,j,k
!	     endif	
             call apex_mall ( &
     	   	gdlat,gdlon,alt,hr,   &		          !Inputs
     	   	b,bhat,bmag,si, 	&		  !Mag Fld
     	   	alon,		      & 		  !Apx Lon
     	   	xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	   	xlatqd,f,f1,f2 ,f3,g1,g2,g3, ist)			 !Qsi-Dpl
              if (ist .ne. 0) stop 'apex ist'
	      fline_r(i,j,isn)%sinI(k) = si          ! sin(I)
	      fline_r(i,j,isn)%D(k)    = d	   	  
	      fline_r(i,j,isn)%F(k)    = f	   	  
!	  
           enddo ! end of height loop
!	 
	   do k=1,fline_p(i,j,isn)%npts 
	     gdlon = fline_p(i,j,isn)%glon(k)  ! geog.longitude should be in [deg]
	     gdlat = fline_p(i,j,isn)%glat(k)  ! geog. latitude should be in [deg]
	   
             hr = h0*1.e-3		      ! reference height convert from [m] to [km]
	     alt =fline_p(i,j,isn)%hgt_pt(k)*1.e-3  ! height convert from [m] to [km]
!	     if(alt < hr ) then
!	   	write(6,*) 'hgt ',hr,alt,i,j,k
!	     endif	
             call apex_mall ( &
     	   	gdlat,gdlon,alt,hr,   &		          !Inputs
     	   	b,bhat,bmag,si, 	&		  !Mag Fld
     	   	alon,		      & 		  !Apx Lon
     	   	xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	   	xlatqd,f,f1,f2 ,f3,g1,g2,g3, ist)			 !Qsi-Dpl
              if (ist .ne. 0) stop 'apex ist'
	      fline_p(i,j,isn)%sinI(k) = si          ! sin(I)
	      fline_p(i,j,isn)%D(k)    = d
	      fline_p(i,j,isn)%F(k)    = f
	      fline_p(i,j,isn)%d1k(k)  = d1(3)  ! k unit upward vector
	      fline_p(i,j,isn)%d2k(k)  = d2(3)  ! k unit upward vector
!		   if(i.eq.70.and.isn.eq.1) write(55,'(2(x,i4),4(x,e15.7))') j,k,alt,d1(3),d1(2),d1(1)
!	      	   	  
!	  
           enddo ! end of height loop
	   
        enddo   ! end latitude loop
!	
!   get values for qd-grid 
!	
	if(isn.eq.1) then ! qd grid: latitude goes from south to northpole
         dlatm =  lat_qd_mp(2)-lat_qd_mp(1)
         do j = 1,nlat_qd-1  ! loop over the midpoints
	   do k=1,nhgt_fix_r
	   
	     if(k.lt.nhgt_fix_r) then
	       qdlat = lat_qd_mp(j)	 ! get quasi dipole latitude midpoint l+/-0.5
	       qdlon = lon_qd_ed(i)	 ! get quasi dipole longitude  edge i-0.5
               hr  = h0*1.e-3		 ! reference height convert from [m] to [km]
	       alt =hgt_qd_mp(k)*1.e-3        ! height convert from [m] to [km] points level-> hgt_fix
               call apex_q2g(qdlat*rtd,qdlon*rtd,alt,gdlat,gdlon,ist)
               if (ist /= 0)  then
	       	  write(6,*) 's2 qdlon/lat/alt= ',qdlon,qdlat,alt
	     	  stop 'apxq2g'
	       end if   
!	   
               call apex_mall ( &
     	   	  gdlat,gdlon,alt,hr,   &		          !Inputs
     	   	  b,bhat,bmag,si, 	&		  !Mag Fld
     	   	  alon,		      & 		  !Apx Lon
     	   	  xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	   	  xlatqd,f,f1,f2 ,f3,g1,g2,g3, ist)			 !Qsi-Dpl
               if (ist .ne. 0) stop 'apex ist'
!       
! M1q  is the area of the qd volume with dlat and dheight (updated 8/3/2015 with July 28, 2015 notes from Art)
!   M1q(i-0.5,l,k) = 2*(r(k)/R)^3 *R^1.5*(r(k+0.5)^0.5-r(k-0.5)^0.5)*dlatqd/F(i-0.5,l,k)
!      with R=Re+h0 and regular spaced in latitude
                fac_r = (hgt_qd_mp(k)+re)/r0  ! r(k)/R
                fac_r = fac_r**3    ! [r(k)/R]**3
                fac_r = 2*fac_r*r0**1.5*(sqrt(hgt_fix_r(k+1)+re)-sqrt(hgt_fix_r(k)+re)) ! 2*[r(k)/R]^3*R^1.5*(r(k+0.5)^0.5-r(k-0.5)^0.5)
                qd(i,j)%M1q(k)= fac_r*dlatm/f
                qd(i,j)%F1(k)  = f
		
	      ! for the i,l,k points I2qd	
	       qdlat = lat_qd_mp(j)	 ! get quasi dipole latitude midpoint l+/-0.5
	       qdlon = lon_qd_mp(i)	 ! get quasi dipole longitude  edge i-0.5
               hr  = h0*1.e-3		 ! reference height convert from [m] to [km]
	       alt =hgt_qd_mp(k)*1.e-3        ! height convert from [m] to [km] points level-> hgt_fix
               call apex_q2g(qdlat*rtd,qdlon*rtd,alt,gdlat,gdlon,ist)
               if (ist /= 0)  then
	       	  write(6,*) 's2 qdlon/lat/alt= ',qdlon,qdlat,alt
	     	  stop 'apxq2g'
	       end if   
!	   
               call apex_mall ( &
     	   	  gdlat,gdlon,alt,hr,   &		          !Inputs
     	   	  b,bhat,bmag,si, 	&		  !Mag Fld
     	   	  alon,		      & 		  !Apx Lon
     	   	  xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	   	  xlatqd,f,f1,f2 ,f3,g1,g2,g3, ist)			 !Qsi-Dpl
               if (ist .ne. 0) stop 'apex ist'
	       
               qd(i,j)%Fqd(k)  = f
               qd(i,j)%f11(k)  = f1(1)
               qd(i,j)%f12(k)  = f1(2)
               qd(i,j)%f21(k)  = f2(1)
               qd(i,j)%f22(k)  = f2(2)
	       
	       if(qdlat.lt.0) then  ! select hemisphere
	          itmp = 1
		  jj = j
	       else
	          itmp = 2
	          jj = nlat_qd-1-j+1 ! one hemisphere on qd midpoints goes from 1 to 80 ; pole to pole is 161
	       endif
		  
	       g13(i,jj,k,itmp) = g1(3)
	       g23(i,jj,k,itmp) = g2(3)
	      
	      endif ! not for k.eq.nhgt_fix_r
!
! M3 for all k levels up to nhgt_fix_r      
	      alt =hgt_qd_ed(k)*1.e-3    ! height convert from [m] to [km] r points half a level below; hgt_qd_ed same as hgt_fix_r
	      qdlat = lat_qd_mp(j)	 ! get quasi dipole latitude midpoint l
	      qdlon = lon_qd_mp(i)	 ! get quasi dipole longitude  midpoints i
!	      
              call apex_q2g(qdlat*rtd,qdlon*rtd,alt,gdlat,gdlon,ist)
              if (ist /= 0)  then
	     	write(6,*) 's2 qdlon/lat/alt= ',qdlon,qdlat,alt
	     	stop 'apxq2g'
	      end if   
              call apex_mall ( &
     	   	gdlat,gdlon,alt,hr,   &		          !Inputs
     	   	b,bhat,bmag,si, 	&		  !Mag Fld
     	   	alon,		      & 		  !Apx Lon
     	   	xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	   	xlatqd,f,f1,f2 ,f3,g1,g2,g3, ist)			 !Qsi-Dpl
              if (ist .ne. 0) stop 'apex ist'
!	      
! M3q is the area of the qd volume with dlat and dlon
! M3q(i,l,k-0.5) = r(k-0.5)^2*dphi_qd/F(i,l,k-0.5)*abs[sin|lam_qd_ed(l-0.5)| - sin|lam_qd_ed(l+0.5)|]
! add an absolute sign since one HP goes from pole to equator, and other goes from equator to pole	      
	      fac_r = (hgt_fix_r(k)+re)**2  ! r(k-0.5)^2
              qd(i,j)%M3q(k)   = fac_r*dlonm/f*abs(sin(abs(lat_qd_ed(j)))-sin(abs(lat_qd_ed(j+1))))
              qd(i,j)%F3(k)    = f
              qd(i,j)%glat(k)  = gdlat
              qd(i,j)%glon(k)  = gdlon
              qd(i,j)%bhat(:,k)= bhat(:)
!	      
            enddo ! end of height loop	   
          enddo   ! end latitude loop
	  
	  endif ! hemisphere isn.eq.1
!	
       enddo   ! end longitude loop
      enddo   ! end of hemipshere loop
!	
     alt = h_leo*1.e-3       ! height convert from [m] to [km] 
     do i=1,nglon
	gdlon = glon(i)  
        do j=1,nglat 
	  gdlat = glat(j)
          call apex_mall ( &
     	     gdlat,gdlon,alt,hr,   &			     !Inputs
     	     b,bhat,bmag,si,	   &		     !Mag Fld
     	     alon,		 &		     !Apx Lon
     	     xlatm,vmp,w,d,be3,sim,d1,d2,d3,e1,e2,e3, & !Mod Apx
     	     xlatqd,f,f1,f2 ,f3,g1,g2,g3, ist)  		    !Qsi-Dpl
          if (ist .ne. 0) stop 'apex ist'
	  bhat_g(i,j,:) = bhat
	  qdlon_g(i,j)  = alon
	  qdlat_g(i,j)  = xlatqd
	  
       enddo   ! end longitude loop
     enddo   ! end longitude loop
!      
      if(debug) write(6,*) 'done apxmall'
!      
      end subroutine apxparm
!--------------------------------------------------------------------------
      subroutine shutdown(msg)
!
! An fatal error has occurred -- shut down the model, including MPI.
!
      implicit none
!
! Args:
      character(len=*) :: msg
!
! Local:
      integer :: ier
      character(len=80) :: errorcode
!
      write(6,"(/,28('>'),' MODEL SHUTDOWN ',28('<'))")
      write(6,"('Shutdown: stop message: ',a)") trim(msg)
!     
      stop 'shutdown'
      end subroutine shutdown
!--------------------------------------------------------------------------------------------
      subroutine cross_product(a,b,c)
!      
!  returns the right-handed vector cross product of two 3-vectors:  C = A x B.    
!      
      implicit none
      
      real, intent(in) :: a(3),b(3)  ! multiplicand & multiplier 3-vector
      real, intent(out) :: c(3)      ! result: 3-vector cross product


      c(1) = a(2)*b(3) - a(3)*b(2)                                                  ! compute cross product components
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)

      end subroutine cross_product
!--------------------------------------------------------------------------------------------
