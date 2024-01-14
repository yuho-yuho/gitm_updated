! EFVM in GITM
! Created: Qingyu Zhu, 08/20/2020
!
! -----------------------------------------------------------------------------

!! GET THE DED1 AND DED2 AT A GITM GRID
! -----------------------------------------------------------------------------
subroutine run_efvm(iBlock)

  use efvm, only: efvm_nmlt, efvm_nmlat, efvm_mlats1, efvm_mlts, &
       efvm_dEd1, efvm_dEd2, efvm_dEd1_sh, efvm_dEd2_sh, &
       efvm_main

  use epm, only: epm_er, epm_disp, epm_scale, epm_scale_sh

  use ModGITM

  implicit none                                                                
                                                                               
  integer, intent(in) :: iBlock                                                
  integer :: iLon, iLat, iAlt                       
  real :: mlat1, mlt1
  real :: dEd1, dEd2
  
  logical :: isNorth = .true.

  dEd1_gitm = 0.
  dEd2_gitm = 0. 

  !! RUN EFVM
  call efvm_main(epm_er,epm_disp,epm_scale,epm_scale_sh)
  
  !! GO THROUGH GITM'S GRIDS
  do iLat = -1,nLats+2 
     do iLon = -1,nLons+2 
        do iAlt = -1,nAlts+2
                                                                               
           ! Get GITM grids                                                    
           mlat1=(MLatitude(iLon, iLat, iAlt, iBlock))                         
           mlt1=MLT(iLon, iLat, iAlt)                                          
                                                                               
           if (mlat1<0) then                                                   
              isNorth = .false.                                                
           else                                                                
              isNorth = .true.                                                 
           end if

           mlat1=abs(mlat1)                                                    
           if (mlt1<0.) mlt1=mlt1+24.                                          
           if (mlt1>24.) mlt1=mlt1-24.

           if (isNorth) then
              call get_grid_value(mlt1,mlat1,efvm_mlts,efvm_mlats1,&   
                   efvm_nmlt,efvm_nmlat,efvm_dEd1,dEd1)
              call get_grid_value(mlt1,mlat1,efvm_mlts,efvm_mlats1,&   
                   efvm_nmlt,efvm_nmlat,efvm_dEd2,dEd2)
           else
              call get_grid_value(mlt1,mlat1,efvm_mlts,efvm_mlats1,&    
                   efvm_nmlt,efvm_nmlat,efvm_dEd1_sh,dEd1)
              call get_grid_value(mlt1,mlat1,efvm_mlts,efvm_mlats1,&    
                   efvm_nmlt,efvm_nmlat,efvm_dEd2_sh,dEd2)
           end if

           dEd1_gitm(iLon,iLat,iAlt) = dEd1 * (1e-3)
           dEd2_gitm(iLon,iLat,iAlt) = dEd2 * (1e-3)

        end do ! Altitude
     end do ! Longitude
  end do ! Latitude

end subroutine run_efvm

!! CONSTRUCT DISTURBED ELECTRIC FIELD IN GITM
! -----------------------------------------------------------------------------
subroutine construct_defield(dir_fac1,dir_fac2,iBlock)

  use ModGITM

  implicit none 

  real, intent(in) :: dir_fac1, dir_fac2
  integer, intent(in) :: iBlock

  integer :: iLon, iLat, iAlt
  real :: dEd1, dEd2 

  do iLat = -1,nLats+2 
     do iLon = -1,nLons+2 
        do iAlt = -1,nAlts+2

           dEd1=dEd1_gitm(iLon,iLat,iAlt)*dir_fac1
           dEd2=dEd2_gitm(iLon,iLat,iAlt)*dir_fac2
           
           dEfield(iLon,iLat,iAlt,:)=dEd1*b0_d1(iLon,iLat,iAlt,:,iBlock)+&
                dEd2*b0_d2(iLon,iLat,iAlt,:,iBlock)

        end do ! Altitude
     end do ! Longitude
  end do ! Latitude
  
end subroutine construct_defield
