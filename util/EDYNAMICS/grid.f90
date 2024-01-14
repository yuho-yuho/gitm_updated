!-------------------------------------------------------------------------------------------
      subroutine gen_highres_grid
!
! grid updated April 2015
! now reference height at k=0.5
!
! grid from Feb. 2015 (see Art's notes 12 Feb. 2015)
! S1 points start at 80 km, closer lat. spacing at low latitudes and in the
!  aurora region
!    
      use params_module, only: nmlat_h,nmlatS2_h,nmlat_T1,nmlat_T2,nmlon,&
        nmlonp1,ylatm,ylonm,ylatm_s,ylonm_s,pi,rho,rho_s,rtd,dtr,re,r0,h0,pi, &
        m2km,km2m,nhgt_fix,nhgt_fix_r,hgt_fix,hgt_fix_r,ha,ha_s
      
      implicit none
!
! Local: 
      integer :: j,jns,k,isn
      real :: b,c,d,e,f,g,h,i,d1,d2,d3,d4,d5,pio2,lam1,lam2,lam3,lam4,lam5, &
        y,y1,y2,y3,y4,y5,lam,th0,lamb,yb,d0,hb, &
        lamc,yc,dc,hs,dhdz,rhob,rhoc,rho1,ymax,ha_loc,rho_loc,hc,h1,r0km,rekm,h0km, &
        cos_avg,dlonm
      real :: apex_height
!      
      parameter (d=6., f=12., h=6.)
      parameter (yb = 5., yc = 14.)
      parameter (d1=25., d2=50., d3=60., d4=72., d5=82.)
      parameter (hs=6, dhdz=6.)    ! [km]
!      
      logical, parameter :: debug=.false.    ! [km]
!
      r0km = r0*m2km
      rekm = re*m2km
      h0km = h0*m2km
      pio2 = pi*0.5              ! pi/2
      hb   = h0km + hs*yb            ! ~ 110.e3 m
      rhob = sqrt(r0km/(hb+rekm))    ! 0.9976829
      lamb = acos(rhob)          ! 0.0680087 [rad]
      hc   = h0km + hs*yc + .5*dhdz*(yc-yb)**2
      rhoc = sqrt(r0km/(hc+rekm))
      lamc = acos(rhoc)
!      
      lam1 = d1*dtr
      lam2 = d2*dtr
      lam3 = d3*dtr
      lam4 = d4*dtr
      lam5 = d5*dtr
!      
      h1 = r0km/cos(lam1)**2 - rekm
      if(debug) write(6,*) 'hb =',hb,'    hc =',hc,'    h1 =',h1
      d0 = lamb*rtd
      dc = lamc*rtd
      if(debug) write(6,*) 'd0 =',d0,'    dc =',dc
      rho1 = cos(lam1)
      if(debug) write(6,*) 'rhob =',rhob,'    rhoc =',rhoc,'    rho1 =',rho1
      c = hs + dhdz*(yc-yb)
      b = (2*r0km*sin(lam1)/(d*rho1**3)+c)/(r0km/rho1**2-rekm-hc)
      if(debug) write(6,*) 'b =',b,'    c =',c
      y1 = yc + alog(b*(h1-hc)/c+1)/b
      if(debug) write(6,*) 'y1 =',y1
      y2 = y1 + d*(lam2-lam1)
      if(debug) write(6,*) 'y2 =',y2
      e = (f-d)/(2.e0*(lam3-lam2))
      y3 = y2 + d*(lam3-lam2) + e*(lam3-lam2)**2
      if(debug) write(6,*) 'y3 =',y3
      y4 = y3 + f*(lam4-lam3)
      if(debug) write(6,*) 'y4 =',y4
      g = (f-h)/(2.e0*(lam5-lam4))
      y5 = y4 + f*(lam5-lam4) - g*(lam5-lam4)**2
      if(debug) write(6,*) 'y5 =',y5
      ymax = y5 + h*(pio2-lam5)
      if(debug) write(6,*) 'ymax =',ymax
      i = (y5+h*(pio2-lam5))/pio2
      if(debug) write(6,10),i
   10 format('i =',e11.4)
      if(debug) write(6,*) 'j    rho         ha     lam(deg)  lam(rad) '
!
      do j=1,nmlat_h   ! goes from equator to pole S1, P points
!	th0 = float(j-1)*pio2/float(nmlat_h-1)
!	y = (j-1)*ymax/float(nmlat_h-1)
! next 2 lines added april 2015 
	th0 = (float(j)-.5)*pio2/(float(nmlat_h)-.5)
	y = (float(j)-.5)*ymax/(float(nmlat_h)-.5)
        ha_loc = h0km + hs*y
        rho_loc = sqrt(r0km/(ha_loc+rekm))
        lam = acos(rho_loc)
	if (y.gt.yb.and.y.le.yc) then 
          ha_loc = h0km + hs*y +.5*dhdz*(y-yb)**2
          rho_loc = sqrt(r0km/(ha_loc+rekm))
          lam = acos(rho_loc)
        endif
	if (y.gt.yc.and.y.le.y1) then
          ha_loc = hc + c*(exp(b*(y-yc))-1)/b
          rho_loc = sqrt(r0km/(ha_loc+rekm))
          lam = acos(rho_loc)
        endif
	if (y.gt.y1.and.y.le.y2) then
       	  lam = (y-y1)/d + lam1
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
	if (y.gt.y2.and.y.le.y3) then 
      	  lam = lam2 + (sqrt(d**2+4.*e*(y-y2))-d)/(2.*e)
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
	if (y.gt.y3.and.y.le.y4) then
      	  lam = (y-y3)/f + lam3
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
	if (y.gt.y4.and.y.le.y5) then 
      	  lam = lam4 + (f-sqrt(f**2-4.*g*(y-y4)))/(2.*g)
          rho_loc = cos(lam)
          ha_loc = r0km/rho_loc**2 - rekm
        endif
	if (y.gt.y5) then
      	  lam = (y-y5)/h + lam5
          rho_loc = cos(lam)
          ha_loc = 9999999.
          if (j.ne.nmlat_h) ha_loc = r0km/rho_loc**2 - rekm
        endif
	
	jns = nmlat_h-j+1     ! jns: nmlat_h (eq) to 1 (pole)
        ylatm(jns,2) =  lam   ! northern hemisphere
        ylatm(jns,1) = -ylatm(jns,2) ! southern hemisphere
        rho(jns,1)   = cos(ylatm(jns,1))  ! cos(ylatm)
        rho(jns,2)   = cos(ylatm(jns,2))  ! cos(ylatm)
! added April 2015	  
! Overwrite numerical inaccuracy in calculating rho for j=nmlat_h / jns = 1: 
        if (jns.eq.1) then
          rho(jns,1)   = 0.
          rho(jns,2)   = 0.
	endif
	ha(jns) = ha_loc*km2m
        if(j.le.nhgt_fix) hgt_fix(j) = ha_loc*km2m
	
	if(debug) write(6,20) jns,rho(jns,2),ha_loc*km2m,90.*ylatm(jns,2)/pio2,ylatm(jns,2)
   20   format(i4,f10.6,2f15.2,f10.4)
!   20   format(i2,f12.8,f11.2,f11.3,f11.5,f11.3,f11.5)
      enddo

! set S2 lat. & rho values 
	if(debug) write(6,*)   " "   
      do j=1,nmlatS2_h
         cos_avg = cos(ylatm(j,2))+cos(ylatm(j+1,2))
         cos_avg = cos_avg*0.5
         ylatm_s(j,2) = acos(cos_avg)
         ylatm_s(j,1) = -ylatm_s(j,2)
         rho_s(j,1) = cos_avg  ! cos(ylatm_s) 
         rho_s(j,2) = cos_avg  ! cos(ylatm_s)
         ha_s(j)    = (r0km/cos_avg**2 - rekm)*km2m
	 if(debug) write(6,'(i4,3(x,f10.6))') j,rho_s(j,2),90.*ylatm_s(j,2)/pio2,ylatm_s(j,2)
	 if(debug) write(44,*) j,ylatm(j,2),ylatm(j+1,2),ylatm_s(j,2)
      enddo 
!      
! set up height levels [m] for r-grid on which Je3 is defined
! these are apex heights of ylatm_s S2 points
      isn=2  ! need to specify only one hemisphere
	if(debug) write(6,*)   " "   
      do k=2,nhgt_fix_r
        j = nmlatS2_h-k+2
        hgt_fix_r(k) = apex_height(rho_s(j,isn))
	if(debug)  write(6,'(i3,1(x,f15.2))') k,hgt_fix_r(k)
      enddo
      k=1
      hgt_fix_r(1) = h0
	if(debug)  write(6,'(i4,3(x,f15.2))') k,hgt_fix_r(k)  
! 
! Magnetic longitudes [rad] 
      dlonm = 2.*pi/float(nmlon) 
      do j=1,nmlon
        ylonm(j) = -pi+float(j-1)*dlonm
      enddo ! i=1,nmlon       
! 
! Magnetic longitudes [rad] for s1 grid which is between the p-points
      do j=1,nmlon
        ylonm_s(j) =ylonm(j)+0.5*dlonm
      enddo ! i=1,nmlon    
!
      end subroutine gen_highres_grid
!-------------------------------------------------------------------------------     
      real function apex_height(rho)
      
      use params_module, only:  re,r0
! calculate apex height of each fieldline
! eq.(3.3) [Richmond 1995]
! ha = R/cos^2lam_m - R_E  with lam_m modified apex latitude
!                               R_E mean Earth radius (code re)
!                               h_R reference height  (code h0)
!                               R = R_E + h_R         (code r0)
      
      implicit none
      real, intent(in):: rho
!      
      apex_height = r0/(rho)**2 -re ! [m]
      
!          
      end function  apex_height
!----------------------------------------------------------------------------
