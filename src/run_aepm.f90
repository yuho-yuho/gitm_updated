! AEPM main in GITM
! Updated from the get_aurora_spectra.f90
! Created: Qingyu Zhu, 05/24/2020
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!! Get the Diff_ef, total_ef, ae from the AEPM at a GITM grid 
! -----------------------------------------------------------------------------
subroutine run_aepm(iBlock)

  use aepm, only: aepm_nmlt, aepm_nmlat, aepm_nchannel, &
       aepm_mlats1, aepm_mlts, aepm_channels, &
       aepm_diff_ef, aepm_diff_ef_sh, &
       aepm_main

  use ModGITM
  use ModInputs, only: ScalingHard, ScalingHardFactor, &
       ScalingSoft, ScalingSoftFactor

  implicit none 

  integer, intent(in) :: iBlock 
  integer :: iLon, iLat, ichannel
  real :: mlat1, mlt1, mlat2, mlt2
  real :: total_ef, total_nf, ae, diff_ef, diff_nf

  logical :: isNorth = .true.

  ElectronEnergyFlux=0.       
  ElectronAverageEnergy=0.
  EleDiffEflux=0.
  EleDiffNflux=0.

  ! Calculate the diff_ef from the AEPM
  call aepm_main

  ! Go through GITM's grids 
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

        do ichannel=1,aepm_nchannel

           if (isNorth) then 
              call get_grid_value(mlt1,mlat1,aepm_mlts,aepm_mlats1,&
                   aepm_nmlt,aepm_nmlat,aepm_diff_ef(:,:,ichannel),diff_ef)
           else
              call get_grid_value(mlt1,mlat1,aepm_mlts,aepm_mlats1,&
                   aepm_nmlt,aepm_nmlat,aepm_diff_ef_sh(:,:,ichannel),diff_ef)
           end if
  

           ! Scaling >500 eV electron differential energy fluxes 
           ! qingyu, 08/30/2020
           if (ScalingHard .and. ichannel<=11) &
                diff_ef = diff_ef * ScalingHardFactor

           ! Scaling <500 eV electron differential energy fluxes
           ! qingyu, 08/30/2020
           if (ScalingSoft .and. (ichannel>11)) &
                diff_ef = diff_ef * ScalingSoftFactor

           ! Get the differential energy/number fluxes
           diff_ef = diff_ef * 1.0e8 * pi
           EleDiffEflux(iLon,iLat,ichannel) = diff_ef
           
           diff_nf = diff_ef/aepm_channels(ichannel)
           EleDiffNflux(iLon,iLat,ichannel) = diff_nf
        end do

        ! Calculate the total energy/number fluxes and average energy
        total_nf=sum(EleDiffEflux(iLon,iLat,1:11))*0.364
        total_ef=dot_product(EleDiffEflux(iLon,iLat,1:11),&
             aepm_channels(1:11))*(1.602*1.0e-12)*0.364

        if (total_nf<=1.0e8) total_nf=1.0e8
        ae=(total_ef/total_nf)/(1.602*1.0e-9) ! keV

        if (total_ef<0.) total_ef=0.
        if (ae<0.01) ae=0.01
        
        ElectronEnergyFlux(iLon,iLat) = total_ef 
        ElectronAverageEnergy(iLon,iLat) = ae


     end do
  end do

end subroutine run_aepm

!!! Use AEPM Diff_nf at 19 energy channels to specify GITM's spectrum 
! -----------------------------------------------------------------------------
subroutine get_aepm_spectra(diff_nf_in,diff_nf_out)

  use aepm, only: aepm_channels, aepm_nchannel
  use ModConstants, only: pi
  use ModSources, only: ED_N_Energies, ED_Energies

  implicit none 

  real, intent(in) :: diff_nf_in(aepm_nchannel) 
  real, intent(out) :: diff_nf_out(ED_N_Energies)

  integer :: i, ii, pos1, pos2
  real :: energy, energy1, energy2, wgt1, wgt2, diff_nf
  real :: channel1(ED_N_Energies), channel2(aepm_nchannel)        
  real :: log_diff_nf_in(aepm_nchannel)        
  real :: min_energy, max_energy


  ! Only use the differential number flux within 30 - 30000 eV
  ! Linearly interpolate in the log space
  ! Note that GITM has more channels than the AEPM

  !! Note that current method may not well interpret the contributions of 
  !! electrons above 30000 eV, a better way to assume a power law for
  !! the high-energy tail

  ! Log 10 of the GITM's channels       
  do i=1,ED_N_Energies                                       
     channel1(i)=log10(ED_Energies(i))
  end do

  do i=1,aepm_nchannel 
     channel2(i)=log10(aepm_channels(i))
     
     if (diff_nf_in(i)>0) then
        log_diff_nf_in(i)=log10(diff_nf_in(i))     
     else          
        log_diff_nf_in(i)=-10.           
     end if  
  end do

  ! Note both channel1 and channel2 are from high to low    
  diff_nf=-10.                         
  min_energy = log10(30.)                                               
  max_energy = log10(30000.) 

  do i=1,ED_N_Energies

     energy = channel1(i)                                                      
                                                                               
     ! > 30000 eV or <30 eV                                                    
     if ((energy>max_energy) .or. (energy<min_energy)) &                       
          diff_nf = -10. 

     ! channel2(1) - 30000 eV                                                  
     if ((energy>=channel2(1)) .and. (energy<=max_energy)) &                   
          diff_nf = log_diff_nf_in(1)                                          
                                                  
     ! 30 eV - channel2(aepm_nchannel)      
     if ((energy<channel2(aepm_nchannel)) .and. (energy>=min_energy)) &      
          diff_nf = log_diff_nf_in(aepm_nchannel)

     ! Other situations                                                        
     if ((energy<channel2(1)) .and. (energy>=channel2(aepm_nchannel))) then  
        
        pos1=-99                                                               
        pos2=-99                                                               
                                                                               
        do ii=1,aepm_nchannel-1                                            
                                                                               
           energy1=channel2(ii)                                                
           energy2=channel2(ii+1) 

           if ((energy<energy1) .and. (energy>=energy2)) then                  
              pos1=ii                                                          
              pos2=ii+1                                                        
              exit                                                            
           end if                                                              
                                                                               
        end do                                                                 
                                                                               
        wgt1=(energy-energy2)/(energy1-energy2)                                
        wgt2=1-abs(wgt1)                                                       
                                                                               
        diff_nf = log_diff_nf_in(pos1)*wgt1 + log_diff_nf_in(pos2)*wgt2

     end if

     diff_nf = 10**diff_nf
     diff_nf_out(i) = diff_nf

  end do

end subroutine get_aepm_spectra
