! Specify GITM's total energy flux and average energy by the 
! self-defined auroral pattern

subroutine get_sdam_aurora(iBlock)

  use ModSDAM
  use ModGITM

  implicit none 

  integer, intent(in) :: iBlock

  integer :: iLon, iLat, iAlt
  real :: mlat1, mlt1, val, val1
  real :: ae

  logical :: isNorth = .true.
  
  ElectronEnergyFlux = 0.
  ElectronAverageEnergy = 0.

  do iLat=-1,nLats+2
     do iLon=-1,nLons+2

        mlat1=(MLatitude(iLon, iLat, nAlts+1, iBlock))
        mlt1=MLT(iLon, iLat, nAlts+1)
        
        if (mlat1<0) then 
           isNorth = .false.                                    
        else                                                      
           isNorth = .true.                                          
        end if

        mlat1=abs(mlat1) 
        if (mlt1<0.) mlt1=mlt1+24.    
        if (mlt1>24.) mlt1=mlt1-24.   

        if (mlat1<50) cycle

        if (isNorth) then                                                      
           call interp_sdam(mlat1,mlt1,sdam_mlats,sdam_mlts,&                  
                sdam_efxin_nh,val)                                             
           call interp_sdam(mlat1,mlt1,sdam_mlats,sdam_mlts,&                  
                sdam_nfxin_nh,val1)
        else
           call interp_sdam(mlat1,mlt1,sdam_mlats,sdam_mlts,&                  
                sdam_efxin_sh,val)                                             
           call interp_sdam(mlat1,mlt1,sdam_mlats,sdam_mlts,&                  
                sdam_nfxin_sh,val1)
        end if

        if (val<=0.) val=0.
        if (val1<1.0e8) val1=1.0e8

        ae=(val/val1)/(1.602*1.0e-9) ! keV
        if (ae<0.01) ae=0.01

        ElectronEnergyFlux(iLon,iLat) = val
        ElectronAverageEnergy(iLon,iLat) = ae

     end do
  end do
end subroutine get_sdam_aurora
