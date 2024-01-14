! This routine is to mapp the field from GITM geo/gmag grids to the grids of
! 3D current model
!
! Created: Qingyu Zhu 07/30/2017, qingyu.zhu@mavs.uta.edu
!
! First, mapping un, vn, sigH, sigP to S1 and S2 grids
! Second, mapping the potential to P grids
! Add in this GITM on 03/02/2020, qingyu  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine gedy_mapfield_s1s2

  ! Map un, vn, sigH, sigP to S1 and S2 grids

  ! Global
  use fieldline_s_module, only: fieldline_s1, fline_s1,&
       fieldline_s2, fline_s2
  
  use params_module, only: nmlon, nmlat_T1, nmlat_h, dtr, h0, &
       nmlat_T2, hgt_fix, rtd, nmlatS2_h, nhgt_fix
  
  use ModGedy, only:  &
       un_gitm, vn_gitm, sigP_gitm, sigH_gitm, &
       lon_gitm, lat_gitm, hgt_gitm, &
       nlon_gitm, nlat_gitm, nhgt_gitm
  
  implicit none

  ! Local
  integer :: ndhgt, ndhgt_top, nhgt_new, istat, &
       lev0, lev1, nlat_gitm2, nmax, ivar_max
  integer :: i, j, k, ii, jj, kk, iStart, isn
  integer, allocatable,dimension(:) :: hgtpos
  integer, allocatable,dimension(:,:,:) :: ip2DS1, ip2DS2, &
       jp2DS1, jp2DS2

  real :: dhgt, dhgt_top, frik, ddlon, ddlat, xloni, &
       friki, frikj, sumP, sumH
  real, allocatable, dimension(:) :: hgt_ext, lat_gitm2
  real, allocatable, dimension(:,:) :: wgt_hgt
  real, allocatable, dimension(:,:,:) :: &
       un_ext, vn_ext, ped_ext, hall_ext, val_ext
  real, allocatable, dimension(:,:,:) :: &
       valhgt, valhgt2, val_s1, val_s2
  real, allocatable, dimension(:,:,:,:) :: &
       wgt2DS1, wgt2DS2

  !!!!!! decay rates of conductivities with height                     
  ! bottom                                                                 
  real, parameter :: rl1_gitm = 5., &
       rl2_gitm = 3.
  ! top                                                                 
  real,parameter :: rtp_gitm = (1./50.+1./150.),&
       rth_gitm = (1./25.+1./150. )

  logical :: found
  
  ! ---------------------------------------------------------------------------
  ! Main
  ! ---------------------------------------------------------------------------

  !!!!!! Extend the size of the height grids
  ! Bottom
  ndhgt= 10

  ! Top (grid scale 16.0)
  dhgt_top = 16.0
  ndhgt_top = 5
  
  ! New size
  nhgt_new = nhgt_gitm + ndhgt + ndhgt_top

  !!!!!! Allocate the new arrays for mapping
  allocate(un_ext(nlon_gitm,nlat_gitm,nhgt_new),stat=istat)
  allocate(vn_ext(nlon_gitm,nlat_gitm,nhgt_new),stat=istat)
  allocate(ped_ext(nlon_gitm,nlat_gitm,nhgt_new),stat=istat)
  allocate(hall_ext(nlon_gitm,nlat_gitm,nhgt_new),stat=istat)
  allocate(val_ext(nlon_gitm,nlat_gitm,nhgt_new),stat=istat)
  allocate(hgt_ext(nhgt_new),stat=istat)

  allocate(hgtpos(nhgt_fix),stat=istat)
  allocate(lat_gitm2(nlat_gitm+2),stat=istat)

  allocate(ip2DS1(nmlon,nmlat_T1,nhgt_fix),stat=istat)
  allocate(jp2DS1(nmlon,nmlat_T1,nhgt_fix),stat=istat)
  allocate(ip2DS2(nmlon,nmlat_T2,nhgt_fix),stat=istat)
  allocate(jp2DS2(nmlon,nmlat_T2,nhgt_fix),stat=istat)

  allocate(wgt_hgt(2,nhgt_fix),stat=istat)
  allocate(wgt2DS1(4,nmlon,nmlat_T1,nhgt_fix),stat=istat)
  allocate(wgt2DS2(4,nmlon,nmlat_T2,nhgt_fix),stat=istat)

  allocate(valhgt(nlon_gitm,nlat_gitm,nhgt_fix),stat=istat)
  allocate(valhgt2(nlon_gitm,nlat_gitm+2,nhgt_fix),stat=istat)
  allocate(val_s1(nmlon,nmlat_T1,nhgt_fix),stat=istat)
  allocate(val_s2(nmlon,nmlat_T2,nhgt_fix),stat=istat)
  
  !!!!!! Fill in the new arrays
  ! Put the original arrays in
  un_ext(:,:,ndhgt+1:nhgt_new-ndhgt_top) = un_gitm(:,:,:)
  vn_ext(:,:,ndhgt+1:nhgt_new-ndhgt_top) = vn_gitm(:,:,:)
  ped_ext(:,:,ndhgt+1:nhgt_new-ndhgt_top) = sigP_gitm(:,:,:)
  hall_ext(:,:,ndhgt+1:nhgt_new-ndhgt_top) = sigH_gitm(:,:,:)
  hgt_ext(ndhgt+1:nhgt_new-ndhgt_top) = hgt_gitm(:)

  !!! Fill the extending values 
  lev0 = ndhgt+1
  lev1 = nhgt_new-ndhgt_top

  ! Bottomside extending values 
  do k= 1, ndhgt
     hgt_ext(k)=h0/1000.0+(k-1)*(hgt_ext(lev0)-h0/1000.0)/&
          (ndhgt*1.0)
     ped_ext(:,:,k)=ped_ext(:,:,lev0)*exp(-(hgt_ext(lev0)-hgt_ext(k))/&
          (rl1_gitm))
     hall_ext(:,:,k)=hall_ext(:,:,lev0)*exp(-(hgt_ext(lev0)-hgt_ext(k))/&
          (rl2_gitm))
     un_ext(:,:,k)=un_ext(:,:,lev0)
     vn_ext(:,:,k)=vn_ext(:,:,lev0)
  enddo

  ! Topside extending values 
  do k=lev1+1,nhgt_new
     hgt_ext(k)=hgt_ext(lev1)+(k-lev1)*dhgt_top
  enddo

  do k=lev1,nhgt_new
     ped_ext(:,:,k)=ped_ext(:,:,lev1-1)*exp(- &
          (hgt_ext(k)-hgt_ext(lev1-1))*rtp_gitm)
     hall_ext(:,:,k)=hall_ext(:,:,lev1-1)*exp(- &
          (hgt_ext(k)-hgt_ext(lev1-1))*rth_gitm)
     un_ext(:,:,k)=un_ext(:,:,lev1-1)
     vn_ext(:,:,k)=vn_ext(:,:,lev1-1)
  enddo

  !!!!!! Calculate the mapping weights
  !!! Vertical first
  iStart = 1
  hgtpos = -1
  do k=1,nhgt_fix
     do i=iStart, nhgt_new-1

        if (hgt_fix(k)/1000.0 .ge. hgt_ext(i) .and. &
             hgt_fix(k)/1000.0 .lt. hgt_ext(i+1)) then
           frik=(hgt_fix(k)/1000.0-hgt_ext(i))/&
                (hgt_ext(i+1)-hgt_ext(i))
           wgt_hgt(1,k)=1.-frik
           wgt_hgt(2,k)=frik
           hgtpos(k)=i
           iStart = i
           exit
        endif
     enddo
  enddo

  !!! Horizontal next (S1, S2 points seperately)
  ip2DS1 = -999
  ip2DS2 = -999
  jp2DS1 = -999
  jp2DS1 = -999

  wgt2DS1 = 0.0
  wgt2DS2 = 0.0

  ddlon = 360.0/nlon_gitm
  ddlat = 180.0/nlat_gitm

  nlat_gitm2 = nlat_gitm +2   ! add pole values              
  lat_gitm2(2:nlat_gitm+1) = lat_gitm
  lat_gitm2(1) = -90.0
  lat_gitm2(nlat_gitm2) = 90.0

  do k=1,nhgt_fix           ! height  
     do i =1, nmlon            ! mlon                    
        do isn=1,2                ! hemisphere

           ! Mapping to S1 points                                         
           do j=1,nmlat_h-k+1         ! S1 points num at each height      
              if(isn .eq. 1) then
                 jj=j
              else
                 jj = nmlat_T1 -j+1
              endif

              ! longitutde                                                
              xloni = (fline_s1(i,j,isn)%glon(k) - lon_gitm(1))/ddlon
              if (xloni<0.) xloni=xloni+dble(nlon_gitm)
              ip2DS1(i,jj,k) = xloni
              friki = xloni- dble(ip2DS1(i,jj,k))
              ip2DS1(i,jj,k)=ip2DS1(i,jj,k)+1

              ! latitude                                                      
              found = .false.
              do kk=1, nlat_gitm2-1
                 if(fline_s1(i,j,isn)%glat(k) < lat_gitm2(kk).or. &
                      fline_s1(i,j,isn)%glat(k) > lat_gitm2(kk+1)) cycle
                 jp2DS1(i,jj,k)=kk
                 frikj = (fline_s1(i,j,isn)%glat(k)-lat_gitm2(kk))/&
                      (lat_gitm2(kk+1)-lat_gitm2(kk))
                 found = .true.
                 exit
              enddo

              wgt2DS1(1,i,jj,k) = (1-friki)*(1-frikj)
              wgt2DS1(2,i,jj,k) =     friki*(1-frikj)
              wgt2DS1(3,i,jj,k) =     friki*frikj
              wgt2DS1(4,i,jj,k) = (1-friki)*frikj
           enddo

           ! Mapping to S2 points                                       
           do j=1,nmlatS2_h-k+1         ! S2 points num at each height        
              if(isn .eq. 1) then
                 jj=j
              else
                 jj = nmlat_T2 -j+1
              endif

              ! longitutde                                              
              xloni = (fline_s2(i,j,isn)%glon(k) - lon_gitm(1))/ddlon
              if (xloni<0.) xloni=xloni+dble(nlon_gitm)
              ip2DS2(i,jj,k) = xloni
              friki = xloni- dble(ip2DS2(i,jj,k))
              ip2DS2(i,jj,k)=ip2DS2(i,jj,k)+1

              ! latitude                                              
              found = .false.
              do kk=1, nlat_gitm2-1
                 if(fline_s2(i,j,isn)%glat(k) < lat_gitm2(kk).or. &
                      fline_s2(i,j,isn)%glat(k) > lat_gitm2(kk+1)) cycle
                 jp2DS2(i,jj,k)=kk
                 frikj = (fline_s2(i,j,isn)%glat(k)-lat_gitm2(kk))/&
                      (lat_gitm2(kk+1)-lat_gitm2(kk))
                 found = .true.
                 exit
              enddo

              wgt2DS2(1,i,jj,k) = (1-friki)*(1-frikj)
              wgt2DS2(2,i,jj,k) =     friki*(1-frikj)
              wgt2DS2(3,i,jj,k) =     friki*frikj
              wgt2DS2(4,i,jj,k) = (1-friki)*frikj
           enddo

        enddo ! hemisphere
     enddo ! mlon
  enddo ! height

  !!!!!! Start mapping different fields 
  
  ivar_max = 4

  do ii=1,ivar_max ! variable

     select case(ii)
     case(1)
        val_ext = ped_ext
     case(2)
        val_ext = hall_ext
     case(3)
        val_ext = un_ext
     case(4)
        val_ext = vn_ext
     case default
        write(*,*) "no such variable"
     end select

     ! 1D vertical mapping (valhgt)                                   
     call map_hgt_gitm(valhgt,val_ext,hgtpos,wgt_hgt,nlon_gitm,&
          nlat_gitm,nhgt_new,nlon_gitm,nlat_gitm,nhgt_fix)

     ! Set the pole values                                       
     ! Avg of all pole -/+1 values                                
     valhgt2(:,2:nlat_gitm+1,:) = valhgt(:,:,:)
     do k =1,nhgt_fix
        valhgt2(:,1,k) = sum(valhgt(:,1,k))/nlon_gitm
        valhgt2(:,nlat_gitm+2,k)=sum(valhgt(:,nlat_gitm,k))/nlon_gitm
     enddo

     ! Horizontal mapping                                               
     ! S1 (val_s1)                                                       
     call map_hor_gitm(val_s1,valhgt2,ip2DS1,jp2DS1,wgt2DS1,nlon_gitm,&
          nlat_gitm2,nmlon,nmlat_T1,nhgt_fix,nmlat_h,nmlat_T1)
     ! S2 (val_s2)                                                        
     call map_hor_gitm(val_s2,valhgt2,ip2DS2,jp2DS2,wgt2DS2,nlon_gitm,&
          nlat_gitm2,nmlon,nmlat_T2,nhgt_fix,nmlatS2_h,nmlat_T2)

     ! copy the values to the fline arrays                                   
     do i = 1, nmlon    ! mlon
        do k=1,nhgt_fix    ! height
           do isn =1,2        ! hemisphere

              ! S1 points come first                                  
              do j=1, nmlat_h-k+1
                 if(isn.eq.1) then
                    jj = j
                 else
                    jj = nmlat_T1 - j + 1
                 endif

                 select case(ii)
                 case(1)
                    fline_s1(i,j,isn)%sigP(k) = val_s1(i,jj,k)
                 case(2)
                    fline_s1(i,j,isn)%sigH(k) = val_s1(i,jj,k)
                 case(3)
                    fline_s1(i,j,isn)%un(k) = val_s1(i,jj,k)
                 case(4)
                    fline_s1(i,j,isn)%vn(k) = val_s1(i,jj,k)
                 case default
                    write(*,*) 'no such variable'
                 end select
              enddo

              ! Next, S2 points                                  
              do j=1, nmlatS2_h-k+1
                 if(isn.eq.1) then
                    jj = j
                 else
                    jj = nmlat_T2 - j + 1
                 endif

                 select case(ii)
                 case(1)
                    fline_s2(i,j,isn)%sigP(k) = val_s2(i,jj,k)
                 case(2)
                    fline_s2(i,j,isn)%sigH(k) = val_s2(i,jj,k)
                 case(3)
                    fline_s2(i,j,isn)%un(k) = val_s2(i,jj,k)
                 case(4)
                    fline_s2(i,j,isn)%vn(k) = val_s2(i,jj,k)
                 case default
                    write(*,*) 'no such variable'
                 end select
              enddo

           enddo  ! hemisphere                                    
        enddo   ! height                                               
     enddo   ! mlon                                                            
  enddo    ! variable 

  !!!!!! Calculate the SigP, SigH
  do i=1,nmlon ! mlon
     do isn=1,2 ! hemisphere

        do j=1,nmlat_h   ! S1 points                                    
           sumP = 0.
           sumH = 0.                                       
           nmax = fline_s1(i,j,isn)%npts
           do k=1,nmax-1
              sumP = sumP + fline_s1(i,j,isn)%sigP(k)*2*&
                   abs(fline_s1(i,j,isn)%Vmp(k+1)-&
                   fline_s1(i,j,isn)%Vmp(k))/&
                   (fline_s1(i,j,isn)%Bmag(k+1)+fline_s1(i,j,isn)%Bmag(k))

              sumH = sumH + fline_s1(i,j,isn)%sigH(k)*2*&
                   abs(fline_s1(i,j,isn)%Vmp(k+1)-&
                   fline_s1(i,j,isn)%Vmp(k))/&
                   (fline_s1(i,j,isn)%Bmag(k+1)+fline_s1(i,j,isn)%Bmag(k))
           enddo
           fline_s1(i,j,isn)%zigP = sumP
           fline_s1(i,j,isn)%zigH = sumH
        enddo

        do j=1,nmlatS2_h   ! S2 points                              
           sumP = 0.
           sumH = 0.                                          
           nmax = fline_s2(i,j,isn)%npts
           do k=1,nmax-1
              sumP = sumP + fline_s2(i,j,isn)%sigP(k)*2*&
                   abs(fline_s2(i,j,isn)%Vmp(k+1)-&
                   fline_s2(i,j,isn)%Vmp(k))/&
                   (fline_s2(i,j,isn)%Bmag(k+1)+fline_s2(i,j,isn)%Bmag(k))

              sumH = sumH + fline_s2(i,j,isn)%sigH(k)*2*&
                   abs(fline_s2(i,j,isn)%Vmp(k+1)-&
                   fline_s2(i,j,isn)%Vmp(k))/&
                   (fline_s2(i,j,isn)%Bmag(k+1)+fline_s2(i,j,isn)%Bmag(k))
           enddo
           fline_s2(i,j,isn)%zigP = sumP
           fline_s2(i,j,isn)%zigH = sumH
        enddo

     enddo ! hemisphere
  enddo ! mlon

end subroutine gedy_mapfield_s1s2

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine gedy_mapfield_p

  ! Map potential to P points

  ! Global
  use ModGedy, only: nmlon_gitm, nmlat_gitm, &
       mlon1_gitm, mlt1_gitm, mlat1_gitm, &
       poten2d_gitm

  use fieldline_p_module, only: fieldline_p, fline_p

  use params_module, only: nmlat_T1, nmlon, ylonm, ylatm, pi, &
       nmlat_h, rtd

  use readhlpoten_module, only: poten_hl

  implicit none

  ! Local
  integer, allocatable, dimension(:) :: ig
  integer :: i, j, k, ii, jj, kk, istart, istat

  real, allocatable, dimension(:,:) :: wt1D, valhalf
  real :: ylat_T1(nmlat_T1), frki
  real :: mlon1_use(nmlon_gitm), mlat1_use(nmlat_gitm) 

  ! ---------------------------------------------------------------------------
  ! Main

  !!! Put ylatm into ylat_T1
  ! Southern hemisphere
  ylat_T1(1:nmlat_h) = ylatm(1:nmlat_h,1)

  ! Northern hemisphere
  do j=1,nmlat_h          ! pole to equator
     jj = nmlat_T1+1 -j    ! equator
     ylat_T1(jj) = ylatm(j,2)
  enddo

  !!! Note that ylonm and ylat_T1 uses rads, 
  ! the best way to reconcile the difference is divide mlon1_gitm
  ! mlat1_gitm by rtd
  mlon1_use = mlon1_gitm/rtd
  mlat1_use = mlat1_gitm/rtd

  !!! Alocate array for mapping
  if (.not.allocated(ig)) allocate(ig(nmlon),stat=istat)
  if (.not.allocated(wt1D)) allocate(wt1D(2,nmlon),stat=istat)     
  if (.not.allocated(valhalf)) allocate(valhalf(nmlon,nmlat_gitm),stat=istat)
  
  !!!!!! Mapping

  !!! First, mapping in longitude
  istart = 1

  ! Since 3D mlon start at -pi, so we need treat it seperately
  i=1
  j=1
  ig(i) = 0
  frki = (mlon1_use(1)-ylonm(1))/(2*pi+mlon1_use(1)-mlon1_use(nmlon_gitm))
  wt1D(1,i) = frki
  wt1D(2,i) = 1.- frki

  do i=2,nmlon 
     do j=istart,nmlon_gitm-1

        if (ylonm(i)>=mlon1_use(j) .and. ylonm(i)<mlon1_use(j+1)) then        
           frki = (mlon1_use(j+1) - ylonm(i))/(mlon1_use(j+1)-mlon1_use(j))
           ig(i) = j
           wt1D(1,i) = frki
           wt1D(2,i) = 1- frki
           istart = j
           exit
        endif

     enddo
  enddo
 
  call map_lon_gitm(valhalf,poten2d_gitm,ig,wt1D,nmlon_gitm,nmlat_gitm, &
       nmlon,nmlat_gitm)

  !!! Second, mapping in latitude
  ! Re-allocate the mapping arrays
  deallocate(ig,stat =istat)
  deallocate(wt1D,stat =istat)

  if (.not.allocated(ig)) allocate(ig(nmlat_T1),stat=istat)
  if (.not.allocated(wt1D)) allocate(wt1D(2,nmlat_T1),stat=istat)

  ! Deal with points on the end
  istart = 1
  i=1
  j=1
  frki = 1.
  ig(i) = j
  wt1D(1,i) = frki
  wt1D(2,i) = 1.- frki
  
  i=nmlat_T1
  j=nmlat_gitm-1
  frki = 0.
  ig(i) = j
  wt1D(1,i) = frki
  wt1D(2,i) = 1.- frki
  
  ! Work on other points
  do i=2,nmlat_T1-1
     do j=istart,nmlat_gitm-1

        if(ylat_T1(i).ge.mlat1_use(j).and.ylat_T1(i).lt.mlat1_use(j+1)) then

           frki = (mlat1_use(j+1)-ylat_T1(i))/(mlat1_use(j+1)-mlat1_use(j))
           ig(i) = j
           wt1D(1,i) = frki
           wt1D(2,i) = 1.- frki
           istart = j
           exit
        endif

     enddo
  enddo

  call map_lat_gitm(poten_hl,valHalf,ig,wt1D,nmlon,nmlat_gitm, &
       nmlon,nmlat_T1)

end subroutine gedy_mapfield_p

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine map_hgt_gitm(f_out,f_in,hgtg,wght,nlon_in,nlat_in,nhgt_in,&
     nlon_out,nlat_out,nhgt_out)
  
  ! Mapping routine in the vertical direction

  ! Arguments
  integer, intent(in) :: nlon_in,nlat_in,nlon_out, &
       nlat_out,nhgt_in,nhgt_out
  integer, intent(in) :: hgtg(nhgt_out)  ! poshgt                     

  real, intent(in) :: f_in(nlon_in,nlat_in,nhgt_in),&
       wght(2,nhgt_out)
  real,intent(out)   :: f_out(nlon_out,nlat_out,nhgt_out)

  ! Local:                                                 
  integer :: k, hgt0, hgt1
  integer :: i,j
  ! ---------------------------------------------------------------------------
  ! Main

  do k=1,nhgt_out
     hgt0=hgtg(k)
     hgt1=hgtg(k)+1
     if (hgt0 == 0) STOP "stop in mag1D_hgt"
     if (hgt1 == nhgt_in+1) STOP "stop in mag1D_hgt"
     f_out(:,:,k) = f_in(:,:,hgt0)*wght(1,k)+f_in(:,:,hgt1)*wght(2,k)
  enddo

end subroutine map_hgt_gitm

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine map_hor_gitm(fmag,fgeo,long,latg,wght,nlon_geo,nlat_geo,&
     nlon_mag,nlat_mag,nhgt,nmlat_half,nmlat_total)

  ! Mapping routine in the horizontal directions

  ! Arguments:                                                  
  integer, intent(in) :: nlon_geo,nlat_geo,nlon_mag,nlat_mag,nhgt,&
       nmlat_half,nmlat_total
  integer, dimension(nlon_mag,nlat_mag,nhgt) :: long, latg

  real,intent(in) :: fgeo(nlon_geo,nlat_geo,nhgt),&
       wght(4,nlon_mag,nlat_mag,nhgt)
  real,intent(out) :: fmag(nlon_mag,nlat_mag,nhgt)

  ! Local:                                                               
  integer :: i, j, k, isn, jj, long0, long1, latg0, latg1

  ! ---------------------------------------------------------------------------
  ! Main

  do k=1,nhgt
     do i=1,nlon_mag
        do j=1,nmlat_half-k+1
           do isn =1,2
              if(isn.eq.1) then
                 jj = j
              else
                 jj = nmlat_total - j + 1
              endif

              if(latg(i,jj,k)+1.gt.nlat_geo) then
                 write(6,*) 'latg > nlat',i,jj,k,latg(i,jj,k),nlat_geo
                 STOP "stop in geo2mag lat.index"
              elseif(latg(i,jj,k).lt.1) then
                 write(6,*) 'latg < 1',i,jj,k,latg(i,jj,k)
                 STOP "stop in geo2mag lat.index"
              endif

              latg0=latg(i,jj,k)
              latg1=latg0+1

                    latg0=latg(i,jj,k)
              latg1=latg0+1

                    long0=long(i,jj,k)
              long1=long0+1
              if (long1 > nlon_geo) long1=1 ! lon:-2.5~2.5                 

              fmag(i,jj,k) = &
                   fgeo(long0,latg0,k)*wght(1,i,jj,k)+&
                   fgeo(long1,latg0,k)*wght(2,i,jj,k)+&
                   fgeo(long1,latg1,k)*wght(3,i,jj,k)+&
                   fgeo(long0,latg1,k)*wght(4,i,jj,k)
              
           enddo
        enddo
     enddo
  enddo

end subroutine map_hor_gitm

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine map_lon_gitm(f_out,f_in,long,wght,nlon_in,nlat_in, &
     nlon_out,nlat_out)

  ! 1D mapping function in longitudinal direction

  ! Arguments
  integer,intent(in) :: nlon_in,nlat_in,nlon_out,nlat_out
  integer,dimension(nlon_out),intent(in) :: long
  real,intent(in) :: f_in(nlon_in,nlat_in), wght(2,nlon_out)
  real,intent(out) :: f_out(nlon_out,nlat_out)

  ! Local
  integer :: i,j,lat,lon0,lon1,lat0,lat1

  ! ---------------------------------------------------------------------------
  ! Main

  do i=1,nlon_out
     lon0 = long(i)
     lon1 = long(i) + 1

     if (lon0 == 0) lon0 = nlon_in         ! no cyclic point               
     if (lon1 == nlon_in+1) lon0 = 1       ! no cyclic point               
     do j=1,nlat_in
        f_out(i,j) =  f_in(lon0,j)*wght(1,i)+f_in(lon1,j)*wght(2,i)
     enddo
  enddo

end subroutine map_lon_gitm

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine map_lat_gitm(f_out,f_in,long,wght,nlon_in,nlat_in, &
     nlon_out,nlat_out)

  ! 1D mapping function in latitude direction

  ! Arguments:
  integer,intent(in) :: nlon_in,nlat_in,nlon_out,nlat_out
  integer,intent(in) :: long(nlat_out)
  real,intent(in)    :: f_in(nlon_in,nlat_in), wght(2,nlat_out)
  real,intent(out)   :: f_out(nlon_out,nlat_out)

  ! Local: 
  integer :: i,j,lat,lon0,lon1,lat0,lat1

  ! ---------------------------------------------------------------------------
  ! Main

  do i=1,nlat_out
     lon0 = long(i)
     lon1 = long(i) + 1
     if (lon0 == 0) STOP "mag1D_lat"
     if (lon1 == nlat_in+1) STOP "mag1D_lat"
     do j=1,nlon_out
        f_out(j,i) =  f_in(j,lon0)*wght(1,i)+f_in(j,lon1)*wght(2,i)
     enddo
  enddo
  
end subroutine map_lat_gitm
