    module lbcurrent_module
!    
! reads in radial current from the lower atmosphere
! upward positive [A/m2]
! mapped to geomagnetic grid, scaled by 1/F and 
! checked that the result integrates to zero over both hemispheres
!    

    use params_module, only: lbcurrent_fname,h0,nmlat_T3,nmlon,nmlat_h, &
        nmlonp1,J3LB
    use fieldline_r_module,only: fieldline_r,fline_r
!     
    implicit none
    
    integer :: nlat_in,nlon_in
    real, allocatable :: JzNo_in(:,:)
    real, allocatable :: CCor_in(:,:)
    real, allocatable :: J_in(:,:)
    real, allocatable :: lon_in(:)
    real, allocatable :: lat_in(:)
!     
    logical, parameter :: debug =.true.	 
!     
    contains
!--------------------------------------------------------------------------------      
    subroutine get_LBcurrent
!    
    use params_module, only: lbcurrent_fname
!      
      implicit none
!     
#include <netcdf.inc>
!    
! these fields will be read in    
!        float Jz_noEPOTEN(time, lat, lon) ;
!                Jz_noEPOTEN:cell_methods = "time: mean" ;
!                Jz_noEPOTEN:long_name = "potential difference" ;
!                Jz_noEPOTEN:units = "V" ;
!        float currentcorr(time, lat, lon) ;
!                currentcorr:units = "A/m2" ;
!                currentcorr:long_name = "GEC source current, corrected for UT" ;
!                currentcorr:cell_methods = "time: mean" ;
!        double lat(lat) ;
!                lat:long_name = "latitude" ;
!                lat:units = "degrees_north" ;
!        double lon(lon) ;
!                lon:long_name = "longitude" ;
!                lon:units = "degrees_east" ; 
!   
! Local 
      integer :: ncid,istat,id,id_nlat,id_nlon,start3(3),count3(3)
!    
! open file
      write(6,*) 'opening nc-file:',lbcurrent_fname
!  open datafile    
      istat = nf_open(lbcurrent_fname,NF_NOWRITE,ncid) 
      if (istat /= NF_NOERR) then
        call check_err(istat,'error w/ open') 
        stop
      endif
!     
! get coordinate dimensions
! nlat_in dimension:
      istat = nf_inq_dimid(ncid,'lat',id_nlat)
      istat = nf_inq_dimlen(ncid,id_nlat,nlat_in)
      if (istat /= NF_NOERR) then
        write(6,"('get_LBcurrent: Error getting nlat dimension from ', &
         'file ',a)") trim(lbcurrent_fname)
        stop
      endif
! nlon_in dimension:
      istat = nf_inq_dimid(ncid,'lon',id_nlon)
      istat = nf_inq_dimlen(ncid,id_nlon,nlon_in)
      if (istat /= NF_NOERR) then
        write(6,"('get_LBcurrent: Error getting nlon dimension from ', &
         'file ',a)") trim(lbcurrent_fname)
        stop
      endif
!      
! allocate arrays
      allocate(JzNo_in(nlon_in,nlat_in),stat =istat)
      if (istat /= 0) STOP "Not enough memory: JzNo_in"
      allocate(CCor_in(nlon_in,nlat_in),stat =istat)
      if (istat /= 0) STOP "Not enough memory: CCor_in"
      allocate(lon_in(nlon_in),stat =istat)
      if (istat /= 0) STOP "Not enough memory: lon_in"
      allocate(lat_in(nlat_in),stat =istat)
      if (istat /= 0) STOP "Not enough memory: lat_in"
      allocate(J_in(nlon_in,nlat_in),stat =istat)
      if (istat /= 0) STOP "Not enough memory: J_in"
!      
! read in coordinates
      istat = nf_inq_varid(ncid,'lat',id)
      if (istat /= NF_NOERR) call check_err(istat,'id lat')
      istat = nf_get_vara_double(ncid,id,1,nlat_in,lat_in)
      if (istat /= NF_NOERR) call check_err(istat,'vara lat')
      
      istat = nf_inq_varid(ncid,'lon',id)
      if (istat /= NF_NOERR) call check_err(istat,'id lon')
      istat = nf_get_vara_double(ncid,id,1,nlon_in,lon_in)
      if (istat /= NF_NOERR) call check_err(istat,'vara lon')
!      
! read_in variables
      start3 =(/1,1,1/)
      count3 =(/nlon_in,nlat_in,1/)
      istat = nf_inq_varid(ncid,'Jz_noEPOTEN',id)
      if (istat /= NF_NOERR) call check_err(istat,'id Jz_noEPOTEN')
      istat = nf_get_vara_double(ncid,id,start3,count3,JzNo_in)
      if (istat /= NF_NOERR) call check_err(istat,'vara Jz_noEPOTEN')
      
      istat = nf_inq_varid(ncid,'currentcorr',id)
      if (istat /= NF_NOERR) call check_err(istat,'id currentcorr')
      istat = nf_get_vara_double(ncid,id,start3,count3,CCor_in)
      if (istat /= NF_NOERR) call check_err(istat,'vara currentcorr')
!      	
      J_in = CCor_in - JzNo_in
!      
! deallocate arrays
      deallocate(JzNo_in,stat =istat)
      if (istat /= 0) STOP "deallocate: JzNo_in"
      deallocate(CCor_in,stat =istat)
      if (istat /= 0) STOP "deallocate: CCor_in"
!      
      end subroutine get_LBcurrent
!--------------------------------------------------------------------------------      
      subroutine map_LBcurrent
!      
      implicit none
! 
! Local:
      integer :: istat,i,j,isn,jj,k,jjg
! 
! set up mapping coefficients
      integer :: istart,iend,idir
      real :: gdlat(nmlon,nmlat_T3),gdlon(nmlon,nmlat_T3), &
        dellon,dellat,xlongi,frki,frkj,Jsum
    integer ::  &
       ig(nmlon,nmlat_T3),   &	  ! geog lon grid containing each geomag point
       jg(nmlon,nmlat_T3)    	  ! geog lat grid containing each geomag point
    real ::  &
       wt(4,nmlon,nmlat_T3), &	  ! weights for geo2mag interpolation
       valmap(nmlon,nmlat_T3)     ! mapped field
!
! Set up parameters for geographic to magnetic interpolation 
! magnetic grid is r grid for I3
! geogr. location for qd points were already calculated in get_apex
!  
!	     fline_r(i,j,isn)%glon(k) = gdlon  ! save geog. longitude [deg]
!	     fline_r(i,j,isn)%glat(k) = gdlat  ! save geog. latitude  [deg]
!  
! put geog. location in global array from hemispheric arrays with isn
      k = 1 ! for the lowest level
      do isn =1,2
        if(isn.eq.1) then ! SH pole to equator
	  istart = 1; iend = nmlat_h; idir = 1
	else  ! NH pole to equator -> change direction
	  istart = nmlat_h; iend = 1; idir = -1
	endif 
!	
        do i=1,nmlon            ! long. loop
	  do j=istart,iend,idir ! lat. loop
	    if(isn.eq.1) jj = j
	    if(isn.eq.2) jj = nmlat_T3 - j + 1
	    gdlon(i,jj) = fline_r(i,j,isn)%glon(k) ! [deg]
	    gdlat(i,jj) = fline_r(i,j,isn)%glat(k) ! [deg]
	  end do
	end do
      end do
!
      dellon = 360./nlon_in  ! assumes regular grid
      dellat = 180./nlat_in  ! assumes regular grid     
      do i=1,nmlon	    ! long. loop
        do j=1,nmlat_T3     ! lat. loop from SPole to NPole
        
        xlongi = (gdlon(i,j) - lon_in(1))/dellon
        if (xlongi < 0.) xlongi = xlongi + float(nlon_in)
        ig(i,j) = xlongi
        frki = xlongi - dble(ig(i,j))
        ig(i,j) = ig(i,j) + 1
        if (ig(i,j) >= nlon_in) ig(i,j) = ig(i,j) - nlon_in
        gdlat = amin1(gdlat,lat_in(nlat_in))
        do jjg=1,nlat_in
          if (gdlat(i,j) > lat_in(jjg)) cycle
          jg(i,j) = jjg - 1
          frkj = (gdlat(i,j) - lat_in(jg(i,j)))/  &
        	 (lat_in(jjg) - lat_in(jg(i,j)))
          jg(i,j) = jg(i,j) + 1
          exit
        enddo
        wt(1,i,j) = (1. - frki)*(1. - frkj)
        wt(2,i,j) =	   frki*(1. - frkj)
        wt(3,i,j) =	   frki*frkj
        wt(4,i,j) = (1. - frki)*frkj
        
        enddo
      enddo 
!          	
! map to geomagnetic grid
      call geo2mag(valmap,J_in,ig,jg,wt,nlon_in,nmlon,nmlon,nmlat_T3)
!          
!  convert valmap back into hemispheric array  
! divide by F at lowest level
!          
      k=1 ! lowest level 
      do isn =1,2
        if(i.eq.1) then ! SH pole to equator
	  istart = 1; iend = nmlat_h; idir = 1
	else  ! NH pole to equator -> change direction
	  istart = nmlat_h; iend = 1; idir = -1
	endif 
!	
        do i=1,nmlon            ! long. loop
	  do j=istart,iend,idir ! lat. loop
	    if(isn.eq.1) jj = j
	    if(isn.eq.2) jj = nmlat_T3 - j + 1
	    J3LB(i,j,isn) = valmap(i,jj) 
	    J3LB(i,j,isn) = J3LB(i,j,isn)/fline_r(i,j,isn)%F(k)  ! divide by scaling factor
	  end do
	end do
      end do     
!          
! check integral_globe = zero
      k=1       ! lowest level
      Jsum = 0.
      do isn =1,2
        do i=1,nmlon 
	  do j=1,nmlat_h
	   Jsum = Jsum + J3LB(i,j,isn)*fline_r(i,j,isn)%M3(k) ! [A/m2*m2]
	  end do
	end do
      end do  
      write(6,*) 'Jsum = ',Jsum   
	      
!          
! deallocate arrays
      deallocate(lon_in,stat =istat)
      if (istat /= 0) STOP "deallocate: lon_in"
      deallocate(lat_in,stat =istat)
      if (istat /= 0) STOP "deallocate: lat_in"		
!      
      end subroutine map_LBcurrent 

!-----------------------------------------------------------------------
      subroutine geo2mag(fmag,fgeo,long,latg,wght,nlonp1_geo,nlat_geo, &
       nlon_mag,nlat_mag)
!
! Transform field fgeo on geographic grid to geomagnetic grid using
!   indices long,latg and weights wght. 
!
! Args:
      integer,intent(in) :: nlonp1_geo,nlon_mag,nlat_mag,nlat_geo
      integer,dimension(nlon_mag,nlat_mag),intent(in) :: long,latg
      real,intent(in) :: fgeo(nlonp1_geo,nlat_geo), &
         wght(4,nlon_mag,nlat_mag)
      real,intent(out) :: fmag(nlon_mag,nlat_mag)
!     integer,intent(in) :: iprint
!
! Local:
      integer :: i,j,lat,lon0,lon1,lat0,lat1
!
      do lat = 1,nlat_mag
        do i=1,nlon_mag
	  lon0 = long(i,lat)
	  lon1 = long(i,lat) + 1
	  lat0 = latg(i,lat)
	  lat1 = latg(i,lat)+1
	  if(lat1 .gt. nlat_geo) lat1 = nlat_geo  ! we are not interested in the
	                                        ! pole anyway
	  if (lon0 == 0) lon0 = nlonp1_geo	! no cyclic point
	  if (lat0 == 0) lat0 = nlat_geo-1	! no cyclic point
          fmag(i,lat) = &
     	   fgeo(lon0,lat0)*wght(1,i,lat)+ &
     	   fgeo(lon1,lat0)*wght(2,i,lat)+ &
     	   fgeo(lon1,lat1)*wght(3,i,lat)+ &
     	   fgeo(lon0,lat1)*wght(4,i,lat)
!
        enddo	
      enddo
!	
      end subroutine geo2mag
!--------------------------------------------------------------------------
      subroutine check_err(status,name)
!      
      implicit none
! 
#include <netcdf.inc>
! 
      integer,intent(in) :: status
      character(len=*),intent(in) :: name
!
      write(6,*) 'error in ',name
      write(6,*) NF_STRERROR(status)
      write(6,*) ' '
!
      stop
! 
      end subroutine check_err
!--------------------------------------------------------------------------------      
      end module lbcurrent_module
