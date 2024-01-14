    module fieldline_s_module
!     
      use params_module, only: r0,re,dtr,ylatm_s,hgt_fix,nhgt_fix, &
        ylonm,ylonm_s,cm2m,nmlon,hgt_fix_r,nhgt_fix_r,rho_s,ha_s
!     
         implicit none
! 
! fieldline_s2 are the points on which Je2 will be calculated (these are
!   in latitudinal direction from the p point   
      type fieldline_s2 
         integer :: npts
!	 
	 real :: ha         ! apex height
	 real :: mlat_m     ! modified apex latitude
	 real :: mlon_m     ! modified apex longitude
	 real :: zigP       ! Pedersen Conductance
	 real :: zigH       ! Hall Conductance
	 real :: Ed1        ! Ed1 electric field
	 real :: Ed2        ! Ed2 electric field
!	 
	 real, pointer :: mlat_qd(:)    ! quasi dipole latitude
	 real, pointer :: mlon_qd(:)    ! quasi dipole longitude
	 real, pointer :: hgt_pt(:)     ! height of point
!	 
	 real, pointer :: glat(:)    ! geog. latitude
	 real, pointer :: glon(:)    ! geog. longitude

	 real, pointer :: Vmp(:)     ! magnetic fpotential [TM]
	 real, pointer :: Bmag(:)    ! magnetic field magnitude [T]
	 real, pointer :: sinI(:)    ! local inclination sin I
	 real, pointer :: bo(:,:)    ! magnetic field component up positive
	 real, pointer :: be3(:)     ! Be3
	 real, pointer :: D(:)       !  D
	 real, pointer :: F(:)       !  F
	 real, pointer :: ds(:)      ! distance between points on fieldline
	 real, pointer :: d1(:,:)    ! d1 vector
	 real, pointer :: d2(:,:)    ! d2 vector
	 real, pointer :: d1d2(:)    ! d1 dot d2 vector
	 real, pointer :: d2d2(:)    ! d2 dot d2 vector
	 real, pointer :: d3(:,:)    ! d3 vector
	 real, pointer :: e1g2(:)    ! e1 dot g2 vector
	 real, pointer :: e2g2(:)    ! e2 dot g2 vector
	 real, pointer :: e1k(:)     ! e1 dot k vector
	 real, pointer :: e2k(:)     ! e2 dot k vector
	 real, pointer :: e3(:,:)    ! e3 vector
	 real, pointer :: M2(:)      ! M1
	 real, pointer :: N2p(:)     ! N1p
	 real, pointer :: N2h(:)     ! N1h
	 real, pointer :: Je2D(:)    ! Je2D
	 real, pointer :: Je2Ion(:)  ! Je2_ionosphere
	 real, pointer :: Je1Ion(:)  ! Je1_ionosphere
	 ! diagnostic
	 real, pointer :: I23d_1(:)  ! I2
	 real, pointer :: I23d_2(:)  ! I2
	 real, pointer :: I23d_3(:)  ! I2
	 ! diagnostic
	 real, pointer :: d1d1(:)  ! dot (d1,d1)
	 real, pointer :: e1g1(:)  ! dot (e1,g1)
	 real, pointer :: e2g1(:)  ! dot (e2,g1)
	 real, pointer :: bg1(:)   ! dot (bhat,g1)
	 real, pointer :: bg2(:)   ! dot (bhat,g2)
	 real, pointer :: Jf1Dyn(:)! Jf1(wind driven)
	 real, pointer :: Jf2Dyn(:)! Jf2(wind driven)
	 
	 real, pointer :: sigH(:)    ! Hall conductivity [S/m]
	 real, pointer :: sigP(:)    ! Pedersen conductivity [S/m]
	 real, pointer :: un(:)      ! zonal neutral wind [m/s] (pos. eastward)
	 real, pointer :: vn(:)      ! meridonal neutral wind [m/s] (pos. northward)
	 ! diagnostic
	 real, pointer :: Ne(:)      ! electron density [#/m3]
	 real, pointer :: Tei(:)     ! Te+Ti [K]
	 	 
	 real, pointer :: je2(:)      ! je2 current
	 real, pointer :: I2(:)       ! I2 current meridional current integrated over the lon/hgt surface
	 
	 integer, pointer :: ngh_pts(:,:)  ! neighboring points lat_index
      end type fieldline_s2
! 
! fieldline_s1 are the points on which Je1 will be calculated (these are
!   in longitudinal direction from the p point   
      type fieldline_s1 
         integer :: npts
!	 
	 real :: ha         ! apex height
	 real :: mlat_m     ! modified apex latitude
	 real :: mlon_m     ! modified apex longitude
	 real :: zigP       ! Pedersen Conductance
	 real :: zigH       ! Hall Conductance
	 real :: Ed1        ! Ed1 electric field
	 real :: Ed2        ! Ed2 electric field
	 real :: ve1        ! ExB ion velocity field
	 real :: ve2        ! ExB ion velocity field
!	 
	 real, pointer :: mlat_qd(:)    ! quasi dipole latitude
	 real, pointer :: mlon_qd(:)    ! quasi dipole longitude
	 real, pointer :: hgt_pt(:)     ! height of point
!	 
	 real, pointer :: glat(:)    ! geog. latitude
	 real, pointer :: glon(:)    ! geog. longitude

	 real, pointer :: Vmp(:)     ! magnetic fpotential [TM]
	 real, pointer :: Bmag(:)    ! magnetic field magnitude []
	 real, pointer :: sinI(:)    ! local inclination sin I
	 real, pointer :: bo(:,:)    ! magnetic field component up positive
	 real, pointer :: be3(:)     ! Be3
	 real, pointer :: D(:)       !  D
	 real, pointer :: F(:)       !  D
	 real, pointer :: ds(:)      ! distance between points on fieldline
	 real, pointer :: d1(:,:)    ! d1 vector
	 real, pointer :: d2(:,:)    ! d2 vector
	 real, pointer :: d3(:,:)    ! d3 vector
	 real, pointer :: d1d1(:)    ! d1 dot d1 vector
	 real, pointer :: d1d2(:)    ! d1 dot d2 vector
	 real, pointer :: d2d2(:)    ! d2 dot d2 vector
	 real, pointer :: e1g2(:)    ! e1 dot g2 vector
	 real, pointer :: e2g2(:)    ! e2 dot g2 vector
	 real, pointer :: e1k(:)     ! e1 dot k vector
	 real, pointer :: e2k(:)     ! e2 dot k vector
	 real, pointer :: e3(:,:)    ! e3 vector
	 real, pointer :: M1(:)      ! M1
	 real, pointer :: N1p(:)     ! N1p
	 real, pointer :: N1h(:)     ! N1h
	 real, pointer :: Je1D(:)    ! Je1D
	 real, pointer :: Je1Ion(:)  ! Je1_ionosphere
	 real, pointer :: Je2Ion(:)  ! Je2_ionosphere
	 real, pointer :: M1Je1D(:)  ! M1Je1D addition to stabilize close to 
!                       	the equator with Cowling conductivity	 
	 real, pointer :: sigH(:)    ! Hall conductivity [S/m]
	 real, pointer :: sigP(:)    ! Pedersen conductivity [S/m]
	 real, pointer :: un(:)      ! zonal neutral wind [m/s] (pos. eastward)
	 real, pointer :: vn(:)      ! meridonal neutral wind [m/s] (pos. northward)
!	 
	 ! diagnostic
	 real, pointer :: Ne(:)      ! electron density [#/m3]
	 real, pointer :: Tei(:)     ! Te+Ti [K]
	 
	 real, pointer :: je1(:)      ! je1 current
	 real, pointer :: I1(:)       ! I1 current eastward current integrated over the lat/hgt surface
	 ! diagnostic
	 real, pointer :: I13d_1(:)  ! I1
	 real, pointer :: I13d_2(:)  ! I1
	 real, pointer :: I13d_3(:)  ! I1
	 
	 integer, pointer :: ngh_pts(:,:)  ! neighboring points lat_index
      end type fieldline_s1
!      
     type (fieldline_s1), allocatable :: fline_s1(:,:,:)
     type (fieldline_s2), allocatable :: fline_s2(:,:,:)
     
     real,parameter ::   Je2Ion_eq(nmlon)=0.      ! at S1 points at k=1 and j=nmlat_h
                                     ! otherwise could be interpolated?
				     ! read in or set later
 !     
     contains
!--------------------------------------------------------------     
      subroutine fieldline_s_dim
           
      use params_module, only: nmlat_h,nmlatS2_h,nmlon,rtd
      use fieldline_p_module,only: fieldline_p,fline_p
!      
      integer :: i,j,k,isn,jns
      integer :: npt_fldline          ! function
      real :: lamqd_from_apex_coord   ! function
!      
      allocate(fline_s2(nmlon,nmlatS2_h,2))   ! note one point less than p-fieldlines
      allocate(fline_s1(nmlon,nmlat_h,2))     ! note same points as p-fieldlines
! 
! relationship between P,S1, and S2 point for the same index (i,j)
!  P(i,j) then is really S1(i+0.5,j) and S2(i,j+0.5) with j increasing equatorward
!  coefficient is calculated at P points
!   
      do isn = 1,2  ! loop over hemisphere
        do j=1,nmlat_h    ! loop over latitudes (direction pole to equator) 
	                  ! s1 point are inbetween p-points with respect to longitude, but ylatm is same as p points	  
      	  fline_s1(:,j,isn)%ha     = fline_p(:,j,isn)%ha               ! apex_height from p-grid 
      	  fline_s1(:,j,isn)%npts   = fline_p(:,j,isn)%npts             ! points on fieldline from p-grid 
      	  fline_s1(:,j,isn)%mlat_m = fline_p(:,j,isn)%mlat_m           ! magnetic latitude  from p-grid 
!        
          do i=1,nmlon ! loop over longitude	    
      	    fline_s1(i,j,isn)%mlon_m = ylonm_s(i)  ! 
      	    allocate(fline_s1(i,j,isn)%hgt_pt(fline_s1(i,j,isn)%npts))  ! should be independent of longitude
      	    allocate(fline_s1(i,j,isn)%mlat_qd(fline_s1(i,j,isn)%npts)) ! should be independent of longitude
      	    allocate(fline_s1(i,j,isn)%mlon_qd(fline_s1(i,j,isn)%npts))
      	    allocate(fline_s1(i,j,isn)%glon(fline_s1(i,j,isn)%npts))
      	    allocate(fline_s1(i,j,isn)%glat(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%Vmp(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%Bmag(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%sinI(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%bo(3,fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%be3(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%D(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%F(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%d1d1(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%d1d2(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%d2d2(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%d1(3,fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%d2(3,fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%d3(3,fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%e1g2(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%e2g2(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%e1k(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%e2k(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%e3(3,fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%M1(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%N1p(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%N1h(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%Je1D(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%M1Je1D(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%Je1Ion(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%Je2Ion(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%sigH(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%sigP(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%un(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%vn(fline_s1(i,j,isn)%npts))
	    !
	    ! diagnostic
	    allocate(fline_s1(i,j,isn)%Ne(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%Tei(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%I13d_1(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%I13d_2(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%I13d_3(fline_s1(i,j,isn)%npts))
	    !
      	    allocate(fline_s1(i,j,isn)%ngh_pts(2,fline_s1(i,j,isn)%npts)) ! lat_ind of neighboring point
	    allocate(fline_s1(i,j,isn)%je1(fline_s1(i,j,isn)%npts))
	    allocate(fline_s1(i,j,isn)%I1(fline_s1(i,j,isn)%npts))
!        	    
      	    do k=1,fline_s1(i,j,isn)%npts
!              
      	      fline_s1(i,j,isn)%hgt_pt(k) = hgt_fix(k)  ! [m] assumes ordering goes from bottom of fieldline to top
      	  					        ! fix heights go also from bottome to top
              fline_s1(i,j,isn)%mlon_qd(k) = ylonm_s(i)   ! independent of latitude and height  
              fline_s1(i,j,isn)%mlat_qd(k) = fline_p(i,j,isn)%mlat_qd(k)   ! quasi dipole latitude from p-grid
      	    enddo
!	    
      	  enddo   ! end loop longitudes
        enddo  ! end loop latitude
	
	
        do j=1,nmlatS2_h    ! loop over latitudes (direction pole to equator) 
	                  ! s2 point are inbetween p-points with respect to rho=cos(ylatm),but same mlon as P points
          fline_s2(:,j,isn)%ha     = ha_s(j)                           ! apex_height
      	  fline_s2(:,j,isn)%npts   = npt_fldline(fline_s2(1,j,isn)%ha) ! points on fieldline
      	  fline_s2(:,j,isn)%mlat_m = ylatm_s(j,isn)                    ! magnetic latitude
!        
          do i=1,nmlon ! loop over longitude
      	    fline_s2(i,j,isn)%mlon_m = ylonm(i)  ! 
      	    
      	    allocate(fline_s2(i,j,isn)%hgt_pt(fline_s2(i,j,isn)%npts))  ! should be independent of longitude
      	    allocate(fline_s2(i,j,isn)%mlat_qd(fline_s2(i,j,isn)%npts)) ! should be independent of longitude
      	    allocate(fline_s2(i,j,isn)%mlon_qd(fline_s2(i,j,isn)%npts))
      	    allocate(fline_s2(i,j,isn)%glon(fline_s2(i,j,isn)%npts))
      	    allocate(fline_s2(i,j,isn)%glat(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Vmp(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Bmag(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%sinI(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%bo(3,fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%be3(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%D(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%F(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%d1d2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%d2d2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%d1(3,fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%d2(3,fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%d3(3,fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%e1g2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%e2g2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%e1k(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%e2k(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%e3(3,fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%M2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%N2p(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%N2h(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Je2D(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Je1Ion(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Je2Ion(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%sigH(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%sigP(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%un(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%vn(fline_s2(i,j,isn)%npts))
	    ! diagnostic
	    allocate(fline_s2(i,j,isn)%Ne(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Tei(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%I23d_1(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%I23d_2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%I23d_3(fline_s2(i,j,isn)%npts))
	    ! diagnostic
	    allocate(fline_s2(i,j,isn)%d1d1(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%e1g1(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%e2g1(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%bg1(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%bg2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Jf1Dyn(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%Jf2Dyn(fline_s2(i,j,isn)%npts))
	    
      	    allocate(fline_s2(i,j,isn)%ngh_pts(2,fline_s2(i,j,isn)%npts)) ! lat_ind of neighboring point
	    allocate(fline_s2(i,j,isn)%je2(fline_s2(i,j,isn)%npts))
	    allocate(fline_s2(i,j,isn)%I2(fline_s2(i,j,isn)%npts))
!        
      	    do k=1,fline_s2(i,j,isn)%npts
!              
      	      fline_s2(i,j,isn)%hgt_pt(k)  = hgt_fix(k)  ! [m] assumes ordering goes from bottom of fieldline to top
      	  					        ! fix heights go also from bottome to top
              fline_s2(i,j,isn)%mlon_qd(k) = ylonm(i)   ! independent of latitude and height 
              fline_s2(i,j,isn)%mlat_qd(k) = lamqd_from_apex_coord(fline_s2(i,j,isn)%mlat_m,hgt_fix(k))   ! quasi dipole latitude
!
      	    enddo
	    
      	  enddo   ! end loop longitudes
        enddo  ! end loop latitude
	
      enddo  ! end loop hemisphere
!      
      end subroutine fieldline_s_dim
!--------------------------------------------------------------------------------      
     end module fieldline_s_module
