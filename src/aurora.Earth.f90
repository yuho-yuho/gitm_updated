subroutine aurora(iBlock)

  use ModGITM
  use ModSources
  use ModTime, only : tSimulation
  use ModInputs
  use ModConstants
  use ModUserGITM

  ! qingyu, 05/24/2020
  use aepm, only: aepm_nchannel

  implicit none

  integer, intent(in) :: iBlock

  real :: alat, hpi, ped, hal, av_kev, eflx_ergs, a,b, maxi
  real :: Factor,temp_ED, avee, eflux, p
  integer :: i, j, k, n, iError, iED
  logical :: IsDone, IsTop, HasSomeAurora

  real, dimension(nLons,nLats,nAlts) :: temp, AuroralBulkIonRate, &
       IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate

  logical :: IsFirstTime(nBlocksMax) = .true.

  real :: f1, f2, f3, f4, f5, power
  real :: de1, de2, de3, de4, de5, detotal, h

  ! qingyu, 05/24/2020
  real :: diff_nf(aepm_nchannel), half_ae

  ! qingyu, 08/12/2020
  real :: ScaleHeight, qqq, rrr, eee, qqq1, ttt

  if (IsFirstTime(iBlock)) then
     IsFirstTime(iBlock) = .false.

     if (iBlock == 1 .and. UseIonPrecipitation) then
        call ReadIonHeat(IonIonizationFilename,  .true.) 
        call ReadIonHeat(IonHeatingRateFilename, .false.) 
     endif

  else
     if (floor((tSimulation - dT)/dTAurora) == &
          floor(tSimulation/dTAurora)) return
  endif

  AuroralBulkIonRate               = 0.0
  AuroralHeatingRate(:,:,:,iBlock) = 0.0

  call report("Aurora",1)
  call start_timing("Aurora")

  if (UseIonPrecipitation) call interpolate_ions( &
       nLons, nLats, nAlts, &
       Longitude(1:nLons,iBlock), Latitude(1:nLats,iBlock), &
       Altitude_GB(1:nLons, 1:nLats, 1:nAlts, iBlock),&
       IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate)

  if (iBlock == 1) then
     HemisphericPowerNorth = 0.0
     HemisphericPowerSouth = 0.0
  endif

  do i=1,nLats
     do j=1,nLons

        UserData2d(j,i,1,2:nUserOutputs,iBlock) = 0.0

        eflx_ergs = ElectronEnergyFlux(j,i) !/ (1.0e-7 * 100.0 * 100.0)
        av_kev    = ElectronAverageEnergy(j,i)

        ! qingyu, 05/24/2020
        diff_nf   = EleDiffNFlux(j,i,:)

        !p = 40.0 * eflx_ergs**0.5 * av_kev / (16.0 + av_kev**2)
        !h = 0.45 * av_kev**0.85 * p
        !write(*,*) "eflux, avee: ", eflx_ergs, av_kev
        !write(*,*) "hall, ped : ", h, p
             
        ! For diffuse auroral models

        ED_Flux = 0.0
        HasSomeAurora = .false.

        if (eflx_ergs > 0.1) then

           UserData2d(j,i,1,2,iBlock) = av_kev
           UserData2d(j,i,1,3,iBlock) = eflx_ergs

!           if (avee > 10.0) write(*,*) "avee, eflux : ",av_kev,eflx_ergs, &
!                j,i,MLatitude(j, i, nAlts+1, iBlock), MLT(j, i, nAlts+1)

           HasSomeAurora = .true.
           avee = av_kev * 1000.0        ! keV -> eV
           eflux = eflx_ergs * 6.242e11  ! ergs/cm2/s -> eV/cm2/s

           power = eflux * Element_Charge*100.0*100.0 * & !(eV/cm2/s -> J/m2/s)
                dLatDist_FB(j, i, nAlts, iBlock) * &
                dLonDist_FB(j, i, nAlts, iBlock)

           if (latitude(i,iBlock) < 0.0) then
              HemisphericPowerSouth = HemisphericPowerSouth + power
           else
              HemisphericPowerNorth = HemisphericPowerNorth + power
           endif

           ! I think that this wrong
           ! a= sqrt(27.0/(2.0*3.14159)) * eflux /(avee**2.5)
           ! The eflux/avee gives the number flux, which is what is the code
           ! needs.
!           a = (eflux/avee) * 2*sqrt(1 / (pi*(avee/2)**3))
           !a = (eflux/avee)* (2**0.5) * ((3/(avee*pi))**(3.0/2.0)) * pi

           !do n=1,ED_N_Energies
              ! I think that this is wrong
              !ED_flux(n) = &
           !a*sqrt(ed_energies(n))*exp(-1.5*ed_energies(n)/avee)

              ! Pat Newell says that while the ratio of the total energy flux
              ! to the number flux is 3kT/2, in reality, it is 2kT, since
              ! the distribution is skewed.
!              ED_flux(n) = &
!                   a*sqrt(ed_energies(n))*exp(-2.0*ed_energies(n)/avee)

           !enddo

           ! qingyu, 05/24/2020

           !! Maxwellian 
           half_ae = avee/2.                                                   
           a = eflux/(2*half_ae**3)                                            
           do n=1,ED_N_Energies                                                
              ED_flux(n) = a * Ed_Energies(n) * exp(-Ed_Energies(n)/half_ae)   
           end do

           !! AEPM spectrum 
           if (UseAEPMAurora .and. UseAEPMSpectra) &
                call get_aepm_spectra(diff_nf,ED_flux)

           ED_flux = ED_flux * pi

        endif

        !!! Newell aurora 
        if (UseNewellAurora .and. UseNewellMono .and. &
             ElectronNumberFluxMono(j,i) > 0.0) then

           av_kev = ElectronEnergyFluxMono(j, i) / &
                    ElectronNumberFluxMono(j, i) * 6.242e11 ! eV

           power = ElectronNumberFluxMono(j, i) * &
                Element_Charge * 100.0 * 100.0 * &    ! (eV/cm2/s -> J/m2/s)
                dLatDist_FB(j, i, nAlts, iBlock) * dLonDist_FB(j, i, nAlts, iBlock)

           if (latitude(i,iBlock) < 0.0) then
              HemisphericPowerSouth = HemisphericPowerSouth + power
           else
              HemisphericPowerNorth = HemisphericPowerNorth + power
           endif

           UserData2d(j,i,1,4,iBlock) = av_kev / 1000.0
           UserData2d(j,i,1,5,iBlock) = ElectronEnergyFluxMono(j, i)

           ! Mono-Energetic goes into one bin only!
           do n=2,ED_N_Energies-1
              if (av_kev < ED_energies(n-1) .and. av_kev >= ED_energies(n)) then
                 ED_flux(n) = ED_Flux(n) + &
                      ElectronNumberFluxMono(j, i) / &
                      (ED_Energies(n-1) - ED_Energies(n))
                 HasSomeAurora = .true.
              endif
           enddo

        endif

        if (UseNewellAurora .and. UseNewellWave .and. &
             ElectronNumberFluxWave(j,i) > 0.0) then

           av_kev = ElectronEnergyFluxWave(j, i) / &
                    ElectronNumberFluxWave(j, i) * 6.242e11 ! eV

           power = ElectronNumberFluxWave(j, i) * &
                Element_Charge * 100.0 * 100.0 * &    ! (eV/cm2/s -> J/m2/s)
                dLatDist_FB(j, i, nAlts, iBlock) * dLonDist_FB(j, i, nAlts, iBlock)

           if (latitude(i,iBlock) < 0.0) then
              HemisphericPowerSouth = HemisphericPowerSouth + power
           else
              HemisphericPowerNorth = HemisphericPowerNorth + power
           endif

           UserData2d(j,i,1,6,iBlock) = av_kev / 1000.0
           UserData2d(j,i,1,7,iBlock) = ElectronEnergyFluxWave(j, i)

           ! Waves goes into five bins only!
           k = 0
           do n=3,ED_N_Energies-3
              if (av_kev < ED_energies(n-1) .and. av_kev >= ED_energies(n)) then
                 k = n
              endif
           enddo
           if (k > 3) then 
              f1 = 1.0
              f2 = 1.2
              f3 = 1.3
              f4 = f2
              f5 = f1
              de1 = ED_energies(k-3)-ED_energies(k-2)
              de2 = ED_energies(k-2)-ED_energies(k-1)
              de3 = ED_energies(k-1)-ED_energies(k)  
              de4 = ED_energies(k)  -ED_energies(k+1)
              de5 = ED_energies(k+1)-ED_energies(k+2)
!              detotal = (de1+de2+de3+de4+de5) * (f1+f2+f3+f4+f5) / 5
              detotal = (f1+f2+f3+f4+f5)
              ED_flux(k-2) = ED_Flux(k-2)+f1*ElectronNumberFluxWave(j, i)/detotal/de1
              ED_flux(k-1) = ED_Flux(k-1)+f2*ElectronNumberFluxWave(j, i)/detotal/de2
              ED_flux(k  ) = ED_Flux(k  )+f3*ElectronNumberFluxWave(j, i)/detotal/de3
              ED_flux(k+1) = ED_Flux(k+1)+f4*ElectronNumberFluxWave(j, i)/detotal/de4
              ED_flux(k+2) = ED_Flux(k+2)+f5*ElectronNumberFluxWave(j, i)/detotal/de5
!              ED_flux(k-2) = ED_Flux(k-2) + f1*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k-1) = ED_Flux(k-1) + f2*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k  ) = ED_Flux(k  ) + f3*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k+1) = ED_Flux(k+1) + f4*ElectronNumberFluxWave(j, i) / detotal
!              ED_flux(k+2) = ED_Flux(k+2) + f5*ElectronNumberFluxWave(j, i) / detotal
           endif

        endif

        if (HasSomeAurora) then

           do n=1,ED_N_Energies
              UserData2d(j,i,1,7+n,iBlock) = ED_flux(n)
           enddo

           call R_ELEC_EDEP (ED_Flux, 15, ED_Energies, 3, ED_Ion, 7)
           call R_ELEC_EDEP (ED_Flux, 15, ED_Energies, 3, ED_Heating, 11)

           iED = 1

!           factor = 1.0
           factor = 0.4
!           factor = 1.0

           do k = 1, nAlts

              p = alog(Pressure(j,i,k,iBlock)*factor)

              IsDone = .false.
              IsTop = .false.
              do while (.not.IsDone)
                 if (ED_grid(iED) >= p .and. ED_grid(iED+1) <= p) then
                    IsDone = .true.
                    ED_Interpolation_Index(k) = iED
                    ED_Interpolation_Weight(k) = (ED_grid(iED) - p) /  &
                         (ED_grid(iED) - ED_grid(iED+1))
                 else
                    if (iED == ED_N_Alts-1) then
                       IsDone = .true.
                       IsTop = .true.
                    else
                       iED = iED + 1
                    endif
                 endif
              enddo

              if (.not.IsTop) then
                 n = ED_Interpolation_Index(k)
                 AuroralBulkIonRate(j,i,k) = ED_Ion(n) - &
                      (ED_Ion(n) - ED_Ion(n+1))*ED_Interpolation_Weight(k)
                 AuroralHeatingRate(j,i,k,iBlock) = ED_Heating(n) - &
                      (ED_Heating(n) - ED_Heating(n+1))*ED_Interpolation_Weight(k)
              else

                 ! Decrease after top of model
                 AuroralBulkIonRate(j,i,k) = ED_Ion(ED_N_Alts) * &
                      factor*Pressure(j,i,k,iBlock) / exp(ED_grid(ED_N_Alts)) 
                 AuroralHeatingRate(j,i,k,iBlock) = ED_Heating(ED_N_Alts) * &
                      factor*Pressure(j,i,k,iBlock) / exp(ED_grid(ED_N_Alts))

              endif

           enddo

           ! Use Fang 2010 Ionization profile
           ! qingyu, 08/12/2020

           if (useFang2010) then

              do k=1,nAlts

                 ! Scale Height [m -> km]
                 ScaleHeight = -Temperature(j,i,k,iBlock)*TempUnit(j,i,k) &  
                      * Boltzmanns_Constant/&                                
                      (Gravity_GB(j,i,k,iBlock)*MeanMajorMass(j,i,k))*1e-3

                 ! Rho [kg/m3 -> g/cm3]
                 rrr = rho(j,i,k,iBlock) * 1.e-3
                 
                 ! Ionization rate [cm3/s]
                 qqq = 0.

                 do n=1,ED_N_Energies 
                 
                    eee = ED_Energies(n)/1000. ! Emono (keV)

                    ! When using previous ionization file, an additional     
                    ! factor is multiplied to ED_Flux                       
                    ! Perhaps not necessary for Fang 2010                   
                    ttt = ED_Flux(n)/pi*(ED_Energies(n)**2)*&               
                         (1.602*1.0e-12)*0.364 ! Qmono

                    call sub_elecmono(eee,ttt,rrr,ScaleHeight,qqq1)

                    ! Sum all channels 
                    qqq = qqq + qqq1

                 end do ! ED_N_Energies

                 ! [cm3/s -> m3/s]
                 AuroralBulkIonRate(j,i,k) = qqq * 1e6

              end do ! Altitude

           end if ! Fang 2010

        endif

     enddo
  enddo

  ! From Rees's book:

  temp = 0.92 * NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock) + &
         1.00 * NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock) + &
         0.56 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)

  if (UseIonPrecipitation) then

     IonPrecipIonRateS(:,:,:,iO_3P_,iBlock)  = &
          0.56*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/temp
     IonPrecipIonRateS(:,:,:,iO2_,iBlock) = &
          1.00*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/temp
     IonPrecipIonRateS(:,:,:,iN2_,iBlock) = &
          0.92*IonPrecipitationBulkIonRate*&
          NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/temp

     IonPrecipHeatingRate(:,:,:,iBlock) = IonPrecipitationHeatingRate

  else

     IonPrecipIonRateS(:,:,:,:,iBlock) = 0.0
     IonPrecipHeatingRate(:,:,:,iBlock) = 0.0

  endif

  AuroralIonRateS(:,:,:,iO_3P_,iBlock)  = &
       0.56*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/temp
  AuroralIonRateS(:,:,:,iO2_,iBlock) = &
       1.00*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iO2_,iBlock)/temp
  AuroralIonRateS(:,:,:,iN2_,iBlock) = &
       0.92*AuroralBulkIonRate*&
       NDensityS(1:nLons,1:nLats,1:nAlts,iN2_,iBlock)/temp

  call end_timing("Aurora")

end subroutine aurora




