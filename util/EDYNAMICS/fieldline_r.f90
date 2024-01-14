     module fieldline_r_module
!     
      use params_module, only: r0,re,dtr,ylatm,hgt_fix_r,hgt_fix,nhgt_fix_r, &
        ylonm,cm2m,rho,ha
!     
         implicit none
!    
      type fieldline_r 
         integer :: npts
!	 
	 real :: ha         ! apex height
	 real :: mlat_m     ! modified apex latitude
	 real :: mlon_m     ! modified apex longitude
	 real :: pot        ! electric potential
!	 
	 real, pointer :: mlat_qd(:)    ! quasi dipole latitude
	 real, pointer :: mlon_qd(:)    ! quasi dipole longitude
	 real, pointer :: hgt_pt(:)     ! height of point
!	 
	 real, pointer :: glat(:)    ! geog. latitude
	 real, pointer :: glon(:)    ! geog. longitude
	 real, pointer :: sinI(:)    ! sinI coefficient
	 real, pointer :: D(:)       ! D coefficient
	 real, pointer :: F(:)       ! F factor
	 real, pointer :: M3(:)      ! M3 coefficient
	 real, pointer :: I3(:)      ! I3 current
	 
	 real, pointer :: a3(:)      ! for mapping: A3(l)/A3(nlat_max)
	 real, pointer :: aa3(:)     ! for mapping from mod.apex to quasi dipole grid A3=Sum_pole_j M3(j)
	 real, pointer :: LI3(:)     ! for mapping from mod.apex to quasi dipole grid latitudinal integrated I3
	 
	 real, pointer :: je3(:)     ! je3 current
	 real, pointer :: Jr(:)      ! Jr current (diagnostic)
	 
	 integer, pointer :: ngh_pts(:,:)  ! neighboring points lat_index
      end type fieldline_r
!      
     type (fieldline_r), allocatable :: fline_r(:,:,:)
!     
     contains
!-------------------------------------------------------------- 
!--------------------------------------------------------------      
      subroutine fieldline_r_dim
          
      use params_module, only: nmlat_h,nmlon
!      
      integer :: i,j,k,nlat_k,lat_k(nmlat_h),isn,jns
      integer :: npt_fldline_r         ! function
      real :: lamqd_from_apex_coord   ! function
!      
      allocate(fline_r(nmlon,nmlat_h,2))
!      
      do isn = 1,2  ! loop over hemisphere
        do j=1,nmlat_h  ! loop over latitudes (pole to equator)
      	  fline_r(:,j,isn)%ha     = ha(j)                              ! apex_height on same fieldline as p points
      	  fline_r(:,j,isn)%npts   = npt_fldline_r(fline_r(1,j,isn)%ha) ! points on fieldline
      	  fline_r(:,j,isn)%mlat_m = ylatm(j,isn)                       ! same as magnetic grid
!        
          do i=1,nmlon ! loop over longitude
      	    fline_r(i,j,isn)%mlon_m = ylonm(i)  ! 
      	    
      	    allocate(fline_r(i,j,isn)%hgt_pt(fline_r(i,j,isn)%npts))  ! should be independent of longitude
      	    allocate(fline_r(i,j,isn)%mlat_qd(fline_r(i,j,isn)%npts)) ! should be independent of longitude
      	    allocate(fline_r(i,j,isn)%mlon_qd(fline_r(i,j,isn)%npts))
      	    allocate(fline_r(i,j,isn)%glon(fline_r(i,j,isn)%npts))
      	    allocate(fline_r(i,j,isn)%glat(fline_r(i,j,isn)%npts))
      	    allocate(fline_r(i,j,isn)%ngh_pts(2,fline_r(i,j,isn)%npts)) ! lat_ind of neighboring point
	    allocate(fline_r(i,j,isn)%D(fline_r(i,j,isn)%npts))
	    allocate(fline_r(i,j,isn)%F(fline_r(i,j,isn)%npts))
	    allocate(fline_r(i,j,isn)%sinI(fline_r(i,j,isn)%npts))
	    allocate(fline_r(i,j,isn)%M3(fline_r(i,j,isn)%npts))
	    allocate(fline_r(i,j,isn)%I3(fline_r(i,j,isn)%npts))
	    
	    allocate(fline_r(i,j,isn)%a3(fline_r(i,j,isn)%npts))
	    allocate(fline_r(i,j,isn)%aa3(fline_r(i,j,isn)%npts))
	    allocate(fline_r(i,j,isn)%LI3(fline_r(i,j,isn)%npts))
	    
	    allocate(fline_r(i,j,isn)%je3(fline_r(i,j,isn)%npts))
	    allocate(fline_r(i,j,isn)%Jr(fline_r(i,j,isn)%npts))
!        
      	    do k=1,fline_r(i,j,isn)%npts
              
      	      fline_r(i,j,isn)%hgt_pt(k) = hgt_fix_r(k)  ! [m] assumes ordering goes from bottom of fieldline to top
      	  					       ! fix heights go also from bottome to top
              fline_r(i,j,isn)%mlon_qd(k) = ylonm(i)   ! independent of latitude and height  
              fline_r(i,j,isn)%mlat_qd(k) = lamqd_from_apex_coord(fline_r(i,j,isn)%mlat_m,hgt_fix_r(k))   ! quasi dipole latitude
!
      	    enddo
      	  enddo
        enddo
      enddo
!
!      
      end subroutine fieldline_r_dim
!--------------------------------------------------------------------------------      
     end module fieldline_r_module
