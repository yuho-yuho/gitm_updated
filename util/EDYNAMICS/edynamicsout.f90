! *******************************************************************
! The output module of Edynamics 
! Adapted from the nclplot.F90
! QZ 10/31/2016
! *******************************************************************

module edynamicsout_module

  ! *****************************************************************
  ! Use other modules
  ! *****************************************************************
      
  use params_module, only: &
       nmlat_h,   & ! number of P,S1,R geomagnetic grid latitudes in one hemisphere
       nmlatS2_h, & ! number of S2 geomagnetic grid latitudes in one hemisphere
       nmlat_T1,  & ! total number of geomagnetic grid latitudes S1 points
       nmlat_T2,  & ! total number of geomagnetic grid latitudes S2 points
       nmlat_T3,  & ! total number of geomagnetic grid latitudes R points at lowest level
       nmlon,&
       nlat_qd,nlat_qd_h, &  !  number of quasi dipole latitudes, edge points j-0.5 and half of hemisphere
       ylatm,ylatm_s,ylonm,ylonm_s,nhgt_fix,hgt_fix,&
       nhgt_fix_r,hgt_fix_r,rtd,&
       J3lb       ! Jr at lb at lower boundary at r points lowest level [A/m2]
       
  use qd_module, only: &
       lat_qd_mp, & ! quasi latitude of midpoint of volume l (nlat_qd)
       lon_qd_mp, &      ! quasi longitude of midpoint of volume (nmlon)
       hgt_qd_ed    ! height of quasi dipole grid = r height level (nhgt_fix_r)
      

  ! *****************************************************************
  ! Parameters which may need to be public
  ! *****************************************************************
      
  integer :: npts1_total,npts2_total,npts3_total,nptsp_total
  integer, parameter :: nmlat2 = nmlat_T2

  logical, parameter:: out_hwm = .true., &
       out_cond = .true., &
       out_nm   = .true., &
       out_jd   = .true.
  
  real :: ylatm2_deg(nmlat2),ylonm2_deg(nmlon), &
       ylatm1_deg(nmlat_T1),ylonm1_deg(nmlon), &
       latqd_deg(nlat_qd-1)
  real, allocatable :: qdlat1_pts(:),hgt1_pts(:), &
       qdlat2_pts(:),hgt2_pts(:), &
       qdlat3_pts(:),hgt3_pts(:), &
       qdlatp_pts(:),hgtp_pts(:)
  real :: qdlat3_lb(nmlat_T3)

  real, allocatable  :: &
       un_fld1(:,:),vn_fld1(:,:),     &
       un_fld2(:,:),vn_fld2(:,:),     &
       sigH_fld1(:,:),sigP_fld1(:,:), &
       sigH_fld2(:,:),sigP_fld2(:,:), &
       sinI_fld1(:,:),sinI_fld2(:,:), &
       N1p_fld1(:,:),N2p_fld2(:,:),   &
       N1h_fld1(:,:),N2h_fld2(:,:),   &
       M1_fld1(:,:),M2_fld2(:,:),M3_fld3(:,:),   &
       je1d_fld1(:,:),je2d_fld2(:,:),S_fldp(:,:), &
       poten_2d(:,:),fac_2d(:,:),zigP1(:,:),zigP2(:,:), &
       fac1_2d(:,:), dfac_2d(:,:), & ! qingyu, 11/25/2020
       poten_nh(:,:), poten_sh(:,:), poten1_2d(:,:), & ! qingyu, 11/26/2020
       Jwd_2d(:,:), &
       zigH1(:,:),zigH2(:,:),   &
       ed1_S1(:,:),ed1_S2(:,:),ve1_S1(:,:), &
       ed2_S1(:,:),ed2_S2(:,:),ve2_S1(:,:), &
       Ex_S1(:,:), Ey_S1(:,:), Ez_S1(:,:), &    !! QZ 07/25
       II1(:,:),II2(:,:),II3(:,:),   &    !! I1, I2, I3 in nclplot
       I1_1(:,:),I1_2(:,:),I1_3(:,:),   &
       I2_1(:,:),I2_2(:,:),I2_3(:,:),   &
       Je1I(:,:),Je2I(:,:),Ne1(:,:),Tei1(:,:),   &
       Ne2(:,:),Tei2(:,:),   &
       Jf1dyn_fld2(:,:),Jf2dyn_fld2(:,:), &
       I3lb_2d(:,:), &
       Jf1_tot(:,:,:),Jf2_tot(:,:,:),Jr_tot(:,:,:),Jeej_tot(:,:,:), &
       Jf1hor_tot(:,:,:),Jf2hor_tot(:,:,:)         
  
   contains
     !==================================================================

     subroutine edynamics_output
 
       use fieldline_s_module,only: fieldline_s1,fline_s1, &
            fieldline_s2,fline_s2
       use fieldline_r_module,only: &
            fieldline_r,fline_r
       use fieldline_p_module,only: &
            fieldline_p,fline_p

       use delB_module,only: delbegrd, delbngrd ,delbugrd, &
            delbe400,delbn400,delbu400,delbscl, &
            delbe70W,delbn70W,delbu70W,psi, &
            qdlon_g,qdlat_g, &
            delbegrd_qd, delbngrd_qd ,delbugrd_qd, &
            delbe400_qd,delbn400_qd,delbu400_qd, &
            delbe0ln_qd,delbn0ln_qd,delbu0ln_qd, &
            delscllH_qd,psi_qd
       
       use qd_module,only: &
            Jf1,&	          ! Jf1(nmlon,nlat_qd_h-1,nhgt_fix,2)
            Jf2,&	          ! Jf2(nmlon,nlat_qd_h-1,nhgt_fix,2)
            Jr,&	          ! Jr(nmlon,nlat_qd_h-1,nhgt_fix,2)
            Jeej,&            ! Jeej(nmlon,nlat_qd_h-1,nhgt_fix,2)
            Jf1hor,&          ! Jf1hor(nmlon,nlat_qd_h-1,nhgt_fix,2)
            Jf2hor            ! Jf2hor(nmlon,nlat_qd_h-1,nhgt_fix,2)     
             
       implicit none
             
       integer :: isn
       integer :: iloop_start,iloop_end,idiff
       integer :: i1,i2,i3,i4,i5,i,j,k
      	
       real :: dlon,dlat	
       
       integer   :: n1,n2,n3,n4,n5,ier
       real, dimension(3) :: npts_d1,npts_d2
       
       !---------------------------------------------------------------

       ! calculate number of total points npts_total on a longitudinal slice
       ! each longitude has the same amount of points

       npts1_total = 0
       npts2_total = 0
       npts3_total = 0
       nptsp_total = 0
       do isn = 1,2
          do j = 1,nmlat_h
             npts1_total =npts1_total + fline_s1(1,j,isn)%npts
             npts3_total =npts3_total + fline_r(1,j,isn)%npts
             nptsp_total =nptsp_total + fline_p(1,j,isn)%npts
          enddo
          do j = 1,nmlatS2_h
             npts2_total =npts2_total + fline_s2(1,j,isn)%npts
          enddo
!	
       enddo

       !---------------------------------------------------------------

       ! put latitude/longitude/height in new matrix
       if (.not. allocated(qdlat1_pts)) then
          allocate(qdlat1_pts(npts1_total))
          allocate(hgt1_pts(npts1_total))
          allocate(qdlat2_pts(npts2_total))
          allocate(hgt2_pts(npts2_total))
          allocate(qdlat3_pts(npts3_total))
          allocate(hgt3_pts(npts3_total))
          allocate(qdlatp_pts(nptsp_total))
          allocate(hgtp_pts(nptsp_total))
       endif
          
       i1=0
       i2=0
       i3=0
       i4=0
       i5=0

       do isn = 1,2
          if(isn.eq.1) then
             iloop_start = 1
             iloop_end   = nmlat_h
             idiff = 1
          elseif(isn.eq.2) then
             iloop_start = nmlat_h
             iloop_end   = 1
             idiff = -1
          endif

          do j=iloop_start,iloop_end,idiff
             do k=1,fline_s1(1,j,isn)%npts   
                ! longitudinal index 1 since all long. are the same

                i1 = i1 +1

                qdlat1_pts(i1)=fline_s1(1,j,isn)%mlat_qd(k)*rtd    
                ! get quasi dipole latitude

                hgt1_pts(i1)  =fline_s1(1,j,isn)%hgt_pt(k)*1.e-3   
                ! height convert from [m] to [km]
             enddo

             do k=1,fline_r(1,j,isn)%npts   
                ! longitudinal index 1 since all long. are the same
                
                i3 = i3 +1
                qdlat3_pts(i3)=fline_r(1,j,isn)%mlat_qd(k)*rtd    
                ! get quasi dipole latitude

                hgt3_pts(i3)  =fline_r(1,j,isn)%hgt_pt(k)*1.e-3   
                ! height convert from [m] to [km]
             enddo

             do k=1,fline_p(1,j,isn)%npts   
                ! longitudinal index 1 since all long. are the same

                i4 = i4 +1
                qdlatp_pts(i4)=fline_p(1,j,isn)%mlat_qd(k)*rtd    
                ! get quasi dipole latitude

                hgtp_pts(i4)  =fline_p(1,j,isn)%hgt_pt(k)*1.e-3   
                ! height convert from [m] to [km]
             enddo

             i5 = i5 +1
             qdlat3_lb(i5)=fline_r(1,j,isn)%mlat_qd(1)*rtd    
             ! get quasi dipole latitude
          enddo
!	
          if(isn.eq.1) then
             iloop_start = 1
             iloop_end   = nmlatS2_h
             idiff = 1
          elseif(isn.eq.2) then
             iloop_start = nmlatS2_h
             iloop_end   = 1
             idiff = -1
          endif

          do j=iloop_start,iloop_end,idiff
             do k=1,fline_s2(1,j,isn)%npts   
                ! longitudinal index 1 since all long. are the same

                i2 = i2 +1
                qdlat2_pts(i2)=fline_s2(1,j,isn)%mlat_qd(k)*rtd   
                ! get quasi dipole latitude

                hgt2_pts(i2)  =fline_s2(1,j,isn)%hgt_pt(k)*1.e-3   
                ! height convert from [m] to [km]
             enddo
          enddo
	
       enddo ! end of hemispheres

       ylatm1_deg(1:nmlat_h) = ylatm(1:nmlat_h,1)
       do j=1,nmlat_h-1
          ylatm1_deg(nmlat_h+j) = ylatm(nmlat_h-j,2)
       end do       
       ylatm1_deg= ylatm1_deg*rtd     !mlat1,potential use
      
       ylatm2_deg(1:nmlatS2_h) = ylatm_s(1:nmlatS2_h,1)
       do j=1,nmlatS2_h
          ylatm2_deg(nmlatS2_h+j) = ylatm_s(nmlatS2_h-j+1,2)
       end do
       ylatm2_deg= ylatm2_deg*rtd     !mlat2
 
       ylonm2_deg= ylonm*rtd   !mlon2, potential use      

       ylonm1_deg= ylonm_s*rtd   !mlon1
       latqd_deg= lat_qd_mp*rtd   !lon_qd
       
       !---------------------------------------------------------------
       ! allocate arrays
       if(out_hwm) then
          if (.not.allocated(un_fld1)) then
             allocate(un_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array un_fld1'
          endif
          if (.not.allocated(vn_fld1)) then
             allocate(vn_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array vn_fld1'
          endif
          if (.not.allocated(un_fld2)) then
             allocate(un_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array un_fld2'
          endif
          if (.not.allocated(vn_fld2)) then
             allocate(vn_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array vn_fld2'
          endif
       endif

       if(out_cond) then
          if (.not.allocated(sigH_fld1)) then
             allocate(sigH_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array sigH_fld1'
          endif
          if (.not.allocated(sigP_fld1)) then
             allocate(sigP_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array sigP_fld1'
          endif
          if (.not.allocated(sigH_fld2)) then
             allocate(sigH_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array sigH_fld2'
          endif
          if (.not.allocated(sigP_fld2)) then
             allocate(sigP_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array sigP_fld2'
          endif
          if (.not.allocated(sinI_fld1)) then
             allocate(sinI_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array sinI_fld1'
          endif
          if (.not.allocated(sinI_fld2)) then
             allocate(sinI_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array sinI_fld2'
          endif
       endif

       if(out_nm) then
          if (.not.allocated(N1p_fld1)) then
             allocate(N1p_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array N1p_fld1'
          endif
          if (.not.allocated(N1h_fld1)) then
             allocate(N1h_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array N1h_fld1'
          endif
          if (.not.allocated(N2h_fld2)) then
             allocate(N2h_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array N2h_fld2'
          endif
          if (.not.allocated(N2p_fld2)) then
             allocate(N2p_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array N2p_fld2'
          endif
          if (.not.allocated(M1_fld1)) then
             allocate(M1_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array M1_fld1'
          endif
          if (.not.allocated(M2_fld2)) then
             allocate(M2_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array M2_fld2'
          endif
          if (.not.allocated(M3_fld3)) then
             allocate(M3_fld3(nmlon,npts3_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array M3_fld3'
          endif
       endif

       if(out_jd) then
          if (.not.allocated(je1d_fld1)) then
             allocate(je1d_fld1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array je1d_fld1'
          endif
          if (.not.allocated(je2d_fld2)) then
             allocate(je2d_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array je2d_fld2'
          endif
          if (.not.allocated(Je1I)) then
             allocate(Je1I(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Je1I'
          endif
          if (.not.allocated(Je2I)) then
             allocate(Je2I(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Je2I'
          endif
          if (.not.allocated(Ne1)) then
             allocate(Ne1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Ne1'
          endif
          if (.not.allocated(Tei1)) then
             allocate(Tei1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Tei1'
          endif
          if (.not.allocated(Ne2)) then
             allocate(Ne2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Ne2'
          endif
          if (.not.allocated(Tei2)) then
             allocate(Tei2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Tei2'
          endif
          if (.not.allocated(S_fldp)) then
             allocate(S_fldp(nmlon,nptsp_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array S_fldp'
          endif
          if (.not.allocated(II1)) then
             allocate(II1(nmlon,npts1_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array I1'
             allocate(I1_1(nmlon,npts1_total),stat=ier)
             allocate(I1_2(nmlon,npts1_total),stat=ier)
             allocate(I1_3(nmlon,npts1_total),stat=ier)
          endif
          if (.not.allocated(II2)) then
             allocate(II2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array I2'
             allocate(I2_1(nmlon,npts2_total),stat=ier)
             allocate(I2_2(nmlon,npts2_total),stat=ier)
             allocate(I2_3(nmlon,npts2_total),stat=ier)
          endif
          if (.not.allocated(II3)) then
             allocate(II3(nmlon,npts3_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array I3'
          endif
          if (.not.allocated(Jf1dyn_fld2)) then
             allocate(Jf1dyn_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Jf1dyn_fld2'
          endif
          if (.not.allocated(Jf2dyn_fld2)) then
             allocate(Jf2dyn_fld2(nmlon,npts2_total),stat=ier)
             if (ier /= 0) write(6,*) 'error allocating array Jf2dyn_fld2'
          endif
          
       endif
       
       if(.not. allocated(poten_2d)) then
          allocate(poten_2d(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array poten_2d'
          allocate(fac_2d(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array fac_2d'
          ! qingyu, 11/25/2020
          allocate(fac1_2d(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array fac1_2d'
          allocate(dfac_2d(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array dfac_2d'

          allocate(poten_nh(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array poten_nh'
          allocate(poten_sh(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array poten_sh'
          allocate(poten1_2d(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array poten1_2d'
          allocate(Jwd_2d(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Jwd_2d'

          allocate(zigP1(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array zigP1'
          allocate(zigP2(nmlon,nmlat_T2),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array zigP2'
          allocate(zigH1(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array zigH1'
          allocate(zigH2(nmlon,nmlat_T2),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array zigH2'
          allocate(ed1_S1(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array ed1_S1'
          allocate(ed2_S1(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array ed2_S1'
          allocate(ve1_S1(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array ve1_S1'
          allocate(ve2_S1(nmlon,nmlat_T1),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array ve2_S1'
          allocate(ed1_S2(nmlon,nmlat_T2),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array ed1_S2'
          allocate(ed2_S2(nmlon,nmlat_T2),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array ed2_S2'
          
          allocate(Ex_S1(nmlon,npts1_total),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Ex_S1'
          allocate(Ey_S1(nmlon,npts1_total),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Ey_S1'
          allocate(Ez_S1(nmlon,npts1_total),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Ez_S1'
          
          allocate(I3lb_2d(nmlon,nmlat_T3),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array I3lb_2d'
          allocate(Jf1_tot(nmlon,nlat_qd-1,nhgt_fix),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Jf1_tot'
          allocate(Jf2_tot(nmlon,nlat_qd-1,nhgt_fix),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Jf2_tot'
          allocate(Jr_tot(nmlon,nlat_qd-1,nhgt_fix),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Jr_tot'
          allocate(Jeej_tot(nmlon,nlat_qd-1,nhgt_fix),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Jeej_tot'
          allocate(Jf1hor_tot(nmlon,nlat_qd-1,nhgt_fix),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Jf1hor_tot'
          allocate(Jf2hor_tot(nmlon,nlat_qd-1,nhgt_fix),stat=ier)
          if (ier /= 0) write(6,*) 'error allocating array Jf2hor_tot'
          
       endif
       !---------------------------------------------------------------

       ! copy the hemispheric array to whole hemisphere
       ! for quasi dipole grid 

       ! southern hemisphere
       isn =1
       Jf1_tot(:,1:nlat_qd_h-1,:)   = Jf1(:,1:nlat_qd_h-1,:,isn)
       Jf2_tot(:,1:nlat_qd_h-1,:)   = Jf2(:,1:nlat_qd_h-1,:,isn)
       Jr_tot(:,1:nlat_qd_h-1,:)    = Jr(:,1:nlat_qd_h-1,:,isn)
       Jeej_tot(:,1:nlat_qd_h-1,:)  = Jeej(:,1:nlat_qd_h-1,:,isn)
       Jf1hor_tot(:,1:nlat_qd_h-1,:)= Jf1hor(:,1:nlat_qd_h-1,:,isn)
       Jf2hor_tot(:,1:nlat_qd_h-1,:)= Jf2hor(:,1:nlat_qd_h-1,:,isn)

       ! northern hemisphere
       isn = 2
       do j=1,nlat_qd_h-1  ! north pole to equator
          n1 = nlat_qd-1-j+1 ! equator to northpole
          Jf1_tot(:,n1,:)   = Jf1(:,j,:,isn)
          Jf2_tot(:,n1,:)   = Jf2(:,j,:,isn)
          Jr_tot(:,n1,:)    = Jr(:,j,:,isn)
          Jeej_tot(:,n1,:)  = Jeej(:,j,:,isn)
          Jf1hor_tot(:,n1,:)= Jf1hor(:,j,:,isn)
          Jf2hor_tot(:,n1,:)= Jf2hor(:,j,:,isn)
       end do
       
       ! for fieldline grid   
       do i=1,nmlon
          n3 = 0
          n2 = 0
          n1 = 0
          n4 = 0
          n5 = 0
          do isn = 1,2
             if(isn.eq.1) then
                iloop_start = 1
                iloop_end   = nmlatS2_h
                idiff = 1
             elseif(isn.eq.2) then
                iloop_start = nmlatS2_h
                iloop_end   = 1
                idiff = -1
             endif

             do j=iloop_start,iloop_end,idiff
                do k=1,fline_s2(i,j,isn)%npts
                   n2 = n2 +1
                   if(out_hwm) then
                      un_fld2(i,n2)=fline_s2(i,j,isn)%un(k)  
                      vn_fld2(i,n2)=fline_s2(i,j,isn)%vn(k)
                   endif
                   if(out_cond) then
                      sigH_fld2(i,n2)=fline_s2(i,j,isn)%sigH(k)  
                      sigP_fld2(i,n2)=fline_s2(i,j,isn)%sigP(k)   
                      sinI_fld2(i,n2)=fline_s2(i,j,isn)%sinI(k) 
                   endif
                   if(out_nm) then
                      N2p_fld2(i,n2)=fline_s2(i,j,isn)%N2p(k)
                      N2h_fld2(i,n2)=fline_s2(i,j,isn)%N2h(k)
                      M2_fld2(i,n2) =fline_s2(i,j,isn)%M2(k)  
                   endif
                   if(out_jd) then
                      je2d_fld2(i,n2)=fline_s2(i,j,isn)%Je2D(k)
                      II2(i,n2)=fline_s2(i,j,isn)%I2(k)
                      I2_1(i,n2)=fline_s2(i,j,isn)%I23d_1(k)
                      I2_2(i,n2)=fline_s2(i,j,isn)%I23d_2(k)
                      I2_3(i,n2)=fline_s2(i,j,isn)%I23d_3(k)
                      Je2I(i,n2)=fline_s2(i,j,isn)%Je2Ion(k)
                      Ne2(i,n2) =fline_s2(i,j,isn)%Ne(k)
                      Tei2(i,n2)=fline_s2(i,j,isn)%Tei(k)
                      Jf1dyn_fld2(i,n2)=fline_s2(i,j,isn)%Jf1Dyn(k)
                      Jf2dyn_fld2(i,n2)=fline_s2(i,j,isn)%Jf2Dyn(k)
                   endif
                enddo

                n5 = j
                if(idiff .eq.-1) n5 =nmlatS2_h+(nmlatS2_h-j+1) 
                zigP2(i,n5) = fline_s2(i,j,isn)%zigP   
                zigH2(i,n5) = fline_s2(i,j,isn)%zigH 
                ed1_S2(i,n5) = fline_s2(i,j,isn)%ed1   
                ed2_S2(i,n5) = fline_s2(i,j,isn)%ed2
             enddo

             if(isn.eq.1) then
                iloop_start = 1
                iloop_end   = nmlat_h
                idiff = 1
             elseif(isn.eq.2) then
                iloop_start = nmlat_h
                iloop_end   = 1
                idiff = -1
             endif

             do j=iloop_start,iloop_end,idiff
                do k=1,fline_s1(i,j,isn)%npts
                   n1 = n1 +1
                   if(out_hwm) then
                      un_fld1(i,n1)=fline_s1(i,j,isn)%un(k)  
                      vn_fld1(i,n1)=fline_s1(i,j,isn)%vn(k)
                   endif
                   if(out_cond) then  
                      sigH_fld1(i,n1)=fline_s1(i,j,isn)%sigH(k)  
                      sigP_fld1(i,n1)=fline_s1(i,j,isn)%sigP(k)  
                      sinI_fld1(i,n1)=fline_s1(i,j,isn)%sinI(k)  
                   endif
                   if(out_nm) then
                      N1p_fld1(i,n1)=fline_s1(i,j,isn)%N1p(k)
                      N1h_fld1(i,n1)=fline_s1(i,j,isn)%N1h(k)
                      M1_fld1(i,n1) =fline_s1(i,j,isn)%M1(k)  
                   endif
                   if(out_jd) then
                      je1d_fld1(i,n1)=fline_s1(i,j,isn)%Je1D(k)
                      II1(i,n1)=fline_s1(i,j,isn)%I1(k)
                      I1_1(i,n1)=fline_s1(i,j,isn)%I13d_1(k)
                      I1_2(i,n1)=fline_s1(i,j,isn)%I13d_2(k)
                      I1_3(i,n1)=fline_s1(i,j,isn)%I13d_3(k)
                      Je1I(i,n1)=fline_s1(i,j,isn)%Je1Ion(k)
                      Ne1(i,n1) =fline_s1(i,j,isn)%Ne(k)
                      Tei1(i,n1)=fline_s1(i,j,isn)%Tei(k)
                   endif
		
                   npts_d1=fline_s1(i,j,isn)%d1(:,k)
                   npts_d2=fline_s1(i,j,isn)%d2(:,k)
                   
                   ! calculate the Ex, Ey, Ez components   !! 07/25
                   
                   Ex_S1(i,n1)=fline_s1(i,j,isn)%ed1*npts_d1(1)+&
                        fline_s1(i,j,isn)%ed2*npts_d2(1)

                   Ey_S1(i,n1)=fline_s1(i,j,isn)%ed1*npts_d1(2)+&	
                        fline_s1(i,j,isn)%ed2*npts_d2(2)		

                   Ez_S1(i,n1)=fline_s1(i,j,isn)%ed1*npts_d1(3)+&
                        fline_s1(i,j,isn)%ed2*npts_d2(3)	        
                enddo

                n5 = j
                if(idiff .eq.-1) n5 =nmlat_h+(nmlat_h-j) 
                zigP1(i,n5) = fline_s1(i,j,isn)%zigP   
                zigH1(i,n5) = fline_s1(i,j,isn)%zigH   
                ed1_S1(i,n5) = fline_s1(i,j,isn)%ed1   
                ed2_S1(i,n5) = fline_s1(i,j,isn)%ed2   
                ve1_S1(i,n5) = fline_s1(i,j,isn)%ve1   
                ve2_S1(i,n5) = fline_s1(i,j,isn)%ve2  

                do k=1,fline_r(i,j,isn)%npts
                   n3 = n3 +1
                   if(out_nm) then
                      M3_fld3(i,n3) =fline_r(i,j,isn)%M3(k)  
                   endif
                   if(out_jd) then
                      II3(i,n3) =fline_r(i,j,isn)%I3(k) 
                   endif
                enddo

                n5 = j
                if(idiff .eq.-1) n5 =nmlat_h+(nmlat_h-j+1) 
                I3lb_2d(i,n5) = J3lb(i,j,isn) 
                do k=1,fline_p(i,j,isn)%npts
                   n4 = n4 +1
                   if(out_jd) then
                      S_fldp(i,n4)=fline_p(i,j,isn)%S(k)
                   endif
                enddo

                n5 = j
                if(idiff .eq.-1) n5 =nmlat_h+(nmlat_h-j) 
                poten_2d(i,n5) = fline_p(i,j,isn)%pot  
                fac_2d(i,n5) = fline_p(i,j,isn)%fac_hl  

                ! qingyu, 11/25/2020
                fac1_2d(i,n5) = fline_p(i,j,isn)%fac_hl2
                dfac_2d(i,n5) = fline_p(i,j,isn)%dfac_hl

                ! qingyu, 11/26/2020
                poten_nh(i,n5) = fline_p(i,j,isn)%pot_nh
                poten_sh(i,n5) = fline_p(i,j,isn)%pot_sh

                ! Combine different components
                poten1_2d(i,n5) = fline_p(i,j,isn)%pot1

                Jwd_2d(i,n5) = (fline_p(i,j,isn)%jwd)/fline_p(i,j,isn)%M3(1) 

             enddo  ! end of latitude loop
             
          enddo ! end of hemisphere loop
       enddo  ! end of longitude loop
      
     end subroutine edynamics_output
     !==================================================================

   end module edynamicsout_module

