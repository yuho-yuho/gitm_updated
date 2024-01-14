     module fieldline_p_module
!     
      use params_module, only: r0,re,dtr,ylatm,hgt_fix,nhgt_fix, &
        ylonm,cm2m,rho,ha
!     
         implicit none
!    
      type fieldline_p 
         integer :: npts
!	 
	 real :: ha         ! apex height
	 real :: mlat_m     ! modified apex latitude
	 real :: mlon_m     ! modified apex longitude
	 
	 real :: pot	   ! electric potential
	 real :: pot_test  ! electric potential am 1/2015 for testing
	 real :: fac_hl	   ! fieldaligned current high latitude
         
         ! qingyu, 11/23/2020 and after
         real :: fac_hl1   ! Read-in FAC
         real :: fac_hl2   ! FAC associated with sysmetric potential
         real :: dfac_hl   ! Residual FAC
         real :: jwd

         ! qingyu, 11/26/2020
         real :: pot_nh    ! Potential correspond to NH residual FAC
         real :: pot_sh    ! Potential correspond to NH residual FAC
         real :: pot1      ! electric potential (symmetric FAC)
!	 
	 real, pointer :: mlat_qd(:)    ! quasi dipole latitude
	 real, pointer :: mlon_qd(:)    ! quasi dipole longitude
	 real, pointer :: hgt_pt(:)     ! height of point
!	 
	 real, pointer :: glat(:)    ! geog. latitude
	 real, pointer :: glon(:)    ! geog. longitude

	 real, pointer :: D(:)       !  D
	 real, pointer :: F(:)       !  F
	 real, pointer :: sinI(:)    !  sinI
	 real, pointer :: d1k(:)     !  d1 dot k vector
	 real, pointer :: d2k(:)     !  d2 dot k vector
	 real, pointer :: M3(:)      !  M3
	 real, pointer :: S(:)       !  S
	 real, pointer :: Jr(:)      !  Jr
	 real, pointer :: I1hor(:)   !  I1horizontal
	 real, pointer :: I2hor(:)   !  I2horizontal
	 
	 integer, pointer :: ngh_pts(:,:)  ! neighboring points lat_index
      end type fieldline_p
!      
     type (fieldline_p), allocatable :: fline_p(:,:,:)
     
     type hgt_fld
         integer :: npts              ! # of latitudes still have points at hgt(k)     
	 integer, pointer :: ilat(:)  ! latitude index
     end type hgt_fld
     type (hgt_fld), allocatable :: hgt_fl(:)	 
!     
     contains
!----------------------------------------------------------------------------- 
      subroutine fieldline_p_dim
           
      use params_module, only: nmlat_h,nmlon
!      
      integer :: i,j,k,nlat_k,lat_k(nmlat_h),isn,ilon,is,jns
      integer :: npt_fldline          ! function
      real :: lamqd_from_apex_coord   ! function
      real :: apex_height             ! function
!      
      allocate(fline_p(nmlon,nmlat_h,2))
!      
      do isn = 1,2  ! loop over hemisphere
        do j=1,nmlat_h  ! loop over latitudes (pole to equator)
      	  fline_p(:,j,isn)%ha     = ha(j)                             ! apex_height
	  fline_p(:,j,isn)%npts   = npt_fldline(fline_p(1,j,isn)%ha)  ! points on fieldline
      	  fline_p(:,j,isn)%mlat_m = ylatm(j,isn)                      ! same as magnetic grid
!        
          do i=1,nmlon ! loop over longitude
      	    fline_p(i,j,isn)%mlon_m = ylonm(i)  ! 
      	    
      	    allocate(fline_p(i,j,isn)%hgt_pt(fline_p(i,j,isn)%npts))  ! should be independent of longitude
      	    allocate(fline_p(i,j,isn)%mlat_qd(fline_p(i,j,isn)%npts)) ! should be independent of longitude
      	    allocate(fline_p(i,j,isn)%mlon_qd(fline_p(i,j,isn)%npts))
      	    allocate(fline_p(i,j,isn)%glon(fline_p(i,j,isn)%npts))
      	    allocate(fline_p(i,j,isn)%glat(fline_p(i,j,isn)%npts))
      	    allocate(fline_p(i,j,isn)%ngh_pts(2,fline_p(i,j,isn)%npts)) ! lat_ind of neighboring point
	    allocate(fline_p(i,j,isn)%D(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%F(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%sinI(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%d1k(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%d2k(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%M3(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%S(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%Jr(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%I1hor(fline_p(i,j,isn)%npts))
	    allocate(fline_p(i,j,isn)%I2hor(fline_p(i,j,isn)%npts))
!	    allocate(fline_p(i,j,isn)%pot(fline_p(i,j,isn)%npts))
!	    allocate(fline_p(i,j,isn)%pot_test(fline_p(i,j,isn)%npts)) ! am 1/2015 for testing
!        
      	    do k=1,fline_p(i,j,isn)%npts
              
      	      fline_p(i,j,isn)%hgt_pt(k) = hgt_fix(k)  ! [m] assumes ordering goes from bottom of fieldline to top
      	  					       ! fix heights go also from bottome to top
              fline_p(i,j,isn)%mlon_qd(k) = ylonm(i)   ! independent of latitude and height  
              fline_p(i,j,isn)%mlat_qd(k) = lamqd_from_apex_coord(fline_p(i,j,isn)%mlat_m,hgt_fix(k))   ! quasi dipole latitude
!
      	    enddo
      	  enddo
        enddo
      enddo
!
! create list with lat at each fixed height
      allocate(hgt_fl(nhgt_fix))
!      
      i=1
      isn = 1 ! southern hemisphere it will be the same in the northern hemisphere
      do k =1, nhgt_fix  ! assumes each longitide is the same (no loop over longitude)
         nlat_k = 0
	 do j=1,nmlat_h  ! latitude loop from pole to equator
	   if(fline_p(i,j,isn)%npts >= k) then           ! check if #of pts on fldline is => height => intersects
	     nlat_k = nlat_k + 1   ! increase number of latitudinal points at that height k
	     lat_k(nlat_k) = j     ! get latitudinal index
	   endif
	 enddo
	 
	 allocate(hgt_fl(k)%ilat(nlat_k))
	 hgt_fl(k)%npts= nlat_k                     ! number of fieldlines intersecting with that height k
	 hgt_fl(k)%ilat(1:nlat_k)= lat_k(1:nlat_k)  ! latitudinal index of fldline intersecting with that height k
	 !
	 ! now use the list of latitudfes at each height to set the neighboring points for each fieldline point
	 do j=1,hgt_fl(k)%npts ! set neighboring points for fldlne
	    do is = 1,2
	    do ilon = 1,nmlon
	    if(j==1) then 
	       fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(1,k)  = -99
	       fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(2,k)  = hgt_fl(k)%ilat(j+1)
	    elseif(j ==  hgt_fl(k)%npts) then
	       fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(1,k)  = hgt_fl(k)%ilat(j-1)
	       fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(2,k)  = -99
	    else
	       fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(1,k)  = hgt_fl(k)%ilat(j-1)
	       fline_p(ilon,hgt_fl(k)%ilat(j),is)%ngh_pts(2,k)  = hgt_fl(k)%ilat(j+1)
	    endif 
	    enddo ! end lon loop
	    enddo   ! end is loop
	 enddo   
      enddo 			
!      
      end subroutine fieldline_p_dim
!--------------------------------------------------------------------------------      
     end module fieldline_p_module
!--------------------------------------------------------------------------------------------
      integer function npt_fldline(apex_height) 
! calculates number of points along a fieldline
! uses apex_height of that fieldline and fixed height grid
!
      use params_module, only: nhgt_fix,hgt_fix
      real,intent(in) :: apex_height
      integer :: i
!      
      i = 1
      npt_fldline =0
      do while (i <= nhgt_fix.and.hgt_fix(i) <= apex_height)  ! round to the nearest number since hgt_fix
        npt_fldline =npt_fldline + 1                                 ! was slightly larger in the last decimal than 90 km
	i = i + 1
	if(i > nhgt_fix) exit
      end do
!            
      end function  npt_fldline
!--------------------------------------------------------------------------------------------
      integer function npt_fldline_r(apex_height) 
! calculates number of points along a fieldline
! uses apex_height of that fieldline and fixed height grid
!
      use params_module, only: nhgt_fix_r,hgt_fix_r
      real,intent(in) :: apex_height
      integer :: i
!      
      i = 1
      npt_fldline_r =0
      do while (i <= nhgt_fix_r.and.hgt_fix_r(i) <= apex_height)  ! round to the nearest number since hgt_fix
        npt_fldline_r =npt_fldline_r + 1                                 ! was slightly larger in the last decimal than 90 km
	i = i + 1
	if(i > nhgt_fix_r) exit
      end do
!            
      end function  npt_fldline_r
!--------------------------------------------------------------------------------      
      real function lamqd_from_apex_coord(latm,hgt_fix)
! calculate quasi dipole latitude lamq from modified apex latitude/mod.apex latitude
!      h height of point
!      lamm mod. apex latitude of the fieldline
! lamq= +/- acos([(Re+h)/(Re+hr)]^0.5*cos(lam_m)) eq. (6.2) [Richmond, 1995]
      use params_module, only: r0,re
!      
      implicit none
!      
      real,intent(in)     :: latm
      real,intent(in)     :: hgt_fix
      real :: fac,tmp  
      real,parameter :: eps = 1.e-6 ! am 10/2014 had to increase for running on glade system/geyser
!      
      fac = (re+hgt_fix)/r0  ! all units [cm]
      fac = sqrt(fac)
      fac = fac*cos(latm)
      
      if (abs(abs(fac)-1.0).lt.eps) fac = 1.0 ! fac was 1.e-15 larger than 1 and cuased problem with acos
! lamq needs to be same sign as lam_m 
      lamqd_from_apex_coord = sign(acos(fac),latm)
! set minimum value to eps
!     if(abs(lamqd_from_apex_coord) <= eps)     lamqd_from_apex_coord =  sign(eps,latm)  
!      
      end function  lamqd_from_apex_coord
!--------------------------------------------------------------------------------------------
