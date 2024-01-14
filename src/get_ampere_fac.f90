! INTERFACE BETWEEN GITM AND GEDY REGARDING THE AMPERE FAC
! ADAPTED FROM THE EARLIER SWARM VERSION
! CREATED: QINGYU ZHU, 11/24/2020
!
! -----------------------------------------------------------------------------

subroutine get_ampere_fac

  use fieldline_p_module, only: fline_p
  use params_module, only: nmlat_h, nmlon

  use ModAMPERE
  use ModGITM
  use ModConstants, only: pi

  implicit none 

  integer :: i, j, isn
  real :: mlonin, mlatin, mltin, fac_out
  real :: mlats(ampr_nmlat), mlts(ampr_nmlt), fac_in(ampr_nmlat,ampr_nmlt)
  real :: fac_in1(ampr_nmlat,ampr_nmlt)
  
  mlats=ampr_mlats_nh(:,1)
  mlts=ampr_mlts_nh(1,:)

  !if (ampr_debug) write(*,*) mlats
  !if (ampr_debug) write(*,*) mlts

  do isn=1,2

     if (isn==1) then
        fac_in(:,:)=ampr_facin_sh(:,:)
     else
        fac_in(:,:)=ampr_facin_nh(:,:)
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

           call interp_fac(mlatin,mltin,mlats,mlts,fac_in,fac_out)
           if (mlatin<50) fac_out=0.
           fline_p(i,j,isn)%fac_hl1 = fac_out*1e-6 ! (muA/m2 --> A/m2)

        end do
     end do

  end do

end subroutine get_ampere_fac
