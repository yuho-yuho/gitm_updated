! INTERFACE BETWEEN GITM AND GEDY REGARDING THE SDFM FAC
! ADAPTED FROM THE EARLIER SWARM VERSION
! CREATED: QINGYU ZHU, 11/24/2020
!
! -----------------------------------------------------------------------------

subroutine get_sdfm_fac

  use fieldline_p_module, only: fline_p
  use params_module, only: nmlat_h, nmlon

  use ModSDFM
  use ModGITM
  use ModConstants, only: pi

  implicit none 

  integer :: i, j, isn
  real :: mlonin, mlatin, mltin, fac_out
  real :: colats(sdfm_nmlat), mlts(sdfm_nmlt), fac_in(sdfm_nmlt,sdfm_nmlat)
  real :: fac_in1(sdfm_nmlt,sdfm_nmlat)
  
  colats=sdfm_colats(:)
  mlts=sdfm_mlts(:)

  !if (sdfm_debug) write(*,*) mlats
  !if (sdfm_debug) write(*,*) mlts

  do isn=1,2

     if (isn==1) then
        fac_in(:,:)=sdfm_facin_sh(:,:)
     else
        fac_in(:,:)=sdfm_facin_nh(:,:)
     end if

     do i=1,nmlon
        do j=1,nmlat_h

           mlatin=fline_p(i,j,isn)%mlat_m/pi*180.
           mlatin=abs(mlatin)
           mlonin=fline_p(i,j,isn)%mlon_m/pi*180.

           call magloctm( &                                                    
                mlonin, &                                                      
                SubsolarLatitude,   &                                          
                SubsolarLongitude,  &                                          
                MagneticPoleColat, &                                           
                MagneticPoleLon,   &                                           
                mltin)

           if (mltin<0) mltin=mltin+24.
           if (mltin>=24) mltin=mltin-24.

           call interp_sdfm(mlatin,mltin,colats,mlts,fac_in,fac_out)
           if (mlatin<50) fac_out=0.
           fline_p(i,j,isn)%fac_hl1 = fac_out*1e-6 ! (muA/m2 --> A/m2)

        end do
     end do

  end do

end subroutine get_sdfm_fac
