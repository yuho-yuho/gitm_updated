! INTERFACE BETWEEN GITM AND GEDY REGARDING THE OIM AMPERE FAC
! ADAPTED FROM GET_AMPERE_FAC.F90
! CREATED: QINGYU ZHU, 12/01/2020
!
! -----------------------------------------------------------------------------

subroutine get_oim_fac

  use fieldline_p_module, only: fline_p
  use params_module, only: nmlat_h, nmlon

  use ModOIM
  use ModGITM
  use ModConstants, only: pi

  implicit none 

  integer :: i, j, isn
  real :: mlonin, mlatin, mltin, fac_out
  real :: mlats(oim_nmlat), mlts(oim_nmlt), fac_in(oim_nmlt,oim_nmlat)
  real :: fac_in1(oim_nmlt,oim_nmlat)

  mlats=oim_mlats
  mlts=oim_mlts

  !if (ampr_debug) write(*,*) mlats
  !if (ampr_debug) write(*,*) mlts

  do isn=1,2

     if (isn==1) then
        fac_in(:,:)=oim_facin_sh(:,:)
     else
        fac_in(:,:)=oim_facin_nh(:,:)
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

           call interp_oimfac(mlatin,mltin,mlats,mlts,fac_in,fac_out)

           ! Smoothed FAC
           ! Unit transformation (muA/m2 --> A/m2)
           if (mlatin<50) fac_out=0.           
           fline_p(i,j,isn)%fac_hl1 = fac_out*1e-6
           !if (oim_debug .and. (i==1)) write(*,*) fac_out

        end do
     end do

  end do

end subroutine get_oim_fac
