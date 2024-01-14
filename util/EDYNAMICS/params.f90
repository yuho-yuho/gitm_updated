!
      module params_module
!      
      implicit none
!
   integer,parameter :: &
! am 10/31/2014 test higher resolution case
!      nmlat_h  = 141,  &          ! a test see grid.f90 / number of magnetic latitudes in one hemisphere P,S1, and R points
       nmlat_h  = 81,   &          ! a number of magnetic latitudes in one hemisphere P,S1, and R points (new grid Feb2015)
!       nmlat_h  = 200,   &          ! b number of magnetic latitudes in one hemisphere P,S1, and R points (new grid Feb2015)
!
       nmlatS2_h= nmlat_h-1,    &  ! number of magnetic latitudes in one hemisphere S2 points
       nmlat_T1 = (2*nmlat_h-1),&  ! total number of magnetic latitudes P/S1 points- equator value is double
       nmlat_T2 = (2*nmlatS2_h),&  ! total number of magnetic latitudes S2 points
       nmlat_T3 = (2*nmlat_h)  ,&  ! total number of magnetic latitudes R points- no equator value
       nmlon    = 100,          &  ! number of magnetic longitudes P,S1,S2,R points
       nmlonp1  = nmlon+1,      &
       nlonlat = nmlon*nmlat_h, &
       nlat_qd = 121,           &  ! number of quasi dipole latitudes, edge points j-0.5
       nlat_qd_h=(nlat_qd+1)*0.5, &! half of the hemisphere (assumes point at equator) edge points j-0.5
!       
! nglon=73, nglat=91 gives 5 x 2 deg lon-lat grid for magnetic perturbations
       nglon = 73, nglat = 91      ! dimension of geographic grid
! 
! lat/lon arrays      
    real :: ylonm(nmlon),    &    ! magnetic longitudes of p and s2 grid; same for both hemispheres
            ylatm(nmlat_h,2),&    ! magnetic latitudes of p and s1 points; 2 index for hemisphere
            ylatm_s(nmlatS2_h,2),&! magnetic latitudes of s2-points; 2 index for hemisphere
            ylonm_s(nmlon),&      ! magnetic longitudes of s1 points; same for both hemispheres
            rho(nmlat_h,2),&      ! cos of magnetic latitudes of p and s1 points; 2 index for hemisphere
            rho_s(nmlatS2_h,2),&  ! cos of magnetic latitudes of s2-points; 2 index for hemisphere 
            glon(nglon),&         ! geographic longitude grid
            glat(nglat)           ! geographic latitude grid
! apex heights for ylam & ylam_s          
    real :: ha(nmlat_h),    &     ! apex height calculated in grid.f90 for ylatm points; same for both hemispheres
            ha_s(nmlatS2_h)       ! apex height calculated in grid.f90 for ylatm_s points; same for both hemispheres
!	       
! for fixed height grid
! am 10/31/2014 test higher resolution case
      integer, parameter :: nhgt_fix   = 51   ! was 51         ! a for feb2015 grid goes up to 760km
!      integer, parameter :: nhgt_fix   =  145 ! for 200 lat for 80 km !87  for 80 km & 141lat! 81   for 90km         ! b for feb2015 grid goes up to 760km
!     integer, parameter :: nhgt_fix   = 55            ! a test see grid.f90/ looked at apex height of ylatm 
!     integer, parameter :: nhgt_fix   = 96           ! b test see grid.f90 /looked at apex height of ylatm 
!     integer, parameter :: nhgt_fix   = 130          ! c test see grid.f90 /looked at apex height of ylatm 
!     integer, parameter :: nhgt_fix   = 130          ! d test see grid.f90 /looked at apex height of ylatm 
!
     integer, parameter :: nhgt_fix_r = nhgt_fix+1   ! looked at apex height of ylatm_s
     real :: hgt_fix(nhgt_fix)  	  ! array with fixed heights for s&p-grod
     real :: hgt_fix_r(nhgt_fix_r)  	  ! array with fixed heights for r-grid (not the highest point?)

       
! globale constants
      real, parameter ::           &
         re = 6.37122e6,           &   ! earth radius (m)		     
         h0 = 8.0e4, r0 =re+h0,    &   ! use mean earth radius [m]
         pi = 3.14159265358979312, &
	 dtr = pi/180.,            &
	 rtd = 180./pi,            &
         cm2m=0.01,                &
         km2m=1.e3,                &
         m2km=1.e-3,               &
         m2cm=100.,                &
	 val_fill=999999.,         &  ! fill value
	 matm = 1.6605e-24,        &  ! mass per atomic weight [g/mole] = 1/N_A (N_A Avodadro number)
	 grav = 8.7,               &  ! gravity [m/s2] 
	 boltz =  1.38065E-23        ! Boltzman constant  [m2 kg/s2/K]
	
	 
! geophysical conditions & time
     real ::  &
       f107,	 &   ! solar radio flux [sfu]
       ap,	 &   ! ap index
       year,	 &   ! year
       doy,	 &   ! day of year
       ut,       &   ! ut
       swden,swvel, & ! solar wind
       bx,by,bz, &    ! IMF conditions (not used so far)
       ctpoten,hpower ! used in tiegcm
	 
! Boundary conditions
! lower boundary Je2LB defined by lower atmosphere model
    logical, parameter :: use_lbJ = .false.
    real  :: J3LB(nmlon,nmlat_h,2)    ! r-points current from lower atmosphere [A/m2]
    character(len=*),parameter :: lbcurrent_fname = & 
       '/glade/p/work/maute/3Dcurrent/LB_data/jaroslav_data.nc'
! 	 
! Boundary conditions
! lower boundary Je2LB defined by lower atmosphere model
    logical, parameter :: test_pot =.false.
!     
! No wind forcing
    logical, parameter :: no_wind =.false.
!   
! TIEGCM readin 
    logical, parameter :: read_tiegcm =.true.
!    
! Empirical conductivities (note does not have auroral conductivity included
!                           therefore do not use with high latitude potential)   
    logical, parameter :: get_empConduc =.false.
! Empirical Wind 
    logical, parameter :: get_empWind =.false.
! 	 
! Jpg
! Ionospheric current
    logical, parameter :: Jpg     =.false.
    logical, parameter :: Jpg_add =.false.  ! flag can be reomved once Jpg is tested

!! QZ 04/24/16
! gitm
! Use Gitm 
   logical, parameter :: usegitm  = .false. 

! 	 
! stabilization Art page 24-26 from Nov 2014
! approximation where net meridional current is small
! 2/5/2015 there still seem to be a bug in there (do NOT USE)
    logical, parameter :: use_stabil =.false. ! DO NOT CHANGE should be false
    
    contains
!-----------------------------------------------------------------
   subroutine set_condition(istep_end)
   
   integer, intent(inout) :: istep_end
 
   if(istep_end.eq.0) then
     f107= 200.     ! solar radio flux [sfu]
     ap  = 6.       ! ap index for use in IRI or MSIS?
     year= 2002.    ! year
     doy = 267.     ! day of year
     ut  = 12.	    ! ut
     istep_end = 1  ! how many time steps 
                    !was 24
   else
     ut = ut + 1.  ! increase time
     if(ut.ge.24) then
       ut  = ut -24.
       doy = doy+1
     endif
   endif  
    
   end subroutine set_condition
!-----------------------------------------------------------------     
   end module params_module
!-------------------------------------------------------------------
