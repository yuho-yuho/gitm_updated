! Subroutines to Get AMIE values on each GITM grid
! Created: Qingyu Zhu, 02/01/2021
! -----------------------------------------------------------------------------

! Potential 
subroutine get_amie_pot(iBlock)

  use ModAMIE
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock

  integer :: iLon, iLat, iAlt
  real :: mlat1, mlt1, val

  logical :: isNorth = .true.

  Potential = 0.
  
  do iLat = -1,nLats+2                                                         
     do iLon = -1,nLons+2                                                      
        do iAlt = -1,nAlts+2

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

           if (mlat1<50) cycle

           if (isNorth) then
              call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&
                   amie_potin_nh,val)
           else
              call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&
                   amie_potin_sh,val)
           end if

           AMIE_Potential(iLon,iLat,iAlt) = val

        end do
     end do
  end do

end subroutine get_amie_pot

! Potential for Gedy
subroutine gedy_amie_pot(nmlt,nmlat,mltin,mlatin,potout)

  use ModAMIE

  implicit none
  
  integer, intent(in) :: nmlt, nmlat                                          
  real, intent(in) :: mltin(nmlt), mlatin(nmlat)                           
  real, intent(out):: potout(nmlt,nmlat) 

  integer :: imlt, imlat                                                      
  real :: mlat1, mlt1, val
  logical :: isNorth = .true.

  do imlt = 1,nmlt   
     do imlat = 1,nmlat 
        
        mlt1=mltin(imlt)                                                       
        mlat1=mlatin(imlat)                                                    
                                                                               
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
           call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&
                amie_potin_nh,val)
        else
           call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&
                amie_potin_sh,val)
        end if

        potout(imlt,imlat) = val
        
     end do
  end do

end subroutine gedy_amie_pot

subroutine get_amie_aurora(iBlock)

  use ModAMIE
  use ModGITM

  ! qingyu, 02/07/2021
  use aepm, only: aepm_calc_HP
  use ModInputs, only: scaleAMIEeflux

  implicit none

  integer, intent(in) :: iBlock

  integer :: iLon, iLat, iAlt
  real :: mlat1, mlt1, val, val1
  real :: NH_HP, SH_HP, scalefac, fac1, fac2 ! qingyu, 02/07/2021

  logical :: isNorth = .true.

  ElectronEnergyFlux = 0.
  ElectronAverageEnergy = 0.

  fac1=1.
  fac2=1.
  
  ! qingyu, 02/07/2021
  ! Get AEPM NH and SH HP
  if (scaleAMIEeflux) then
     call aepm_calc_HP(NH_HP,SH_HP)
     fac1=NH_HP/amie_hpiin_nh
     fac2=SH_HP/amie_hpiin_sh
  end if

  do iLat = -1,nLats+2                                                         
     do iLon = -1,nLons+2                                                      

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
           call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&
                amie_efxin_nh,val)
           call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&              
                amie_ekvin_nh,val1)
           scalefac=fac1
        else
           call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&
                amie_efxin_sh,val)
           call interp_amie(mlat1,mlt1,amie_mlats,amie_mlts,&              
                amie_ekvin_sh,val1)           
           scalefac=fac2
        end if

        ElectronEnergyFlux(iLon,iLat) = val * scalefac
        ElectronAverageEnergy(iLon,iLat) = val1

     end do
  end do
  
end subroutine get_amie_aurora
