! New chemistry module developed by Cissi Lin
! Added in this GITM by qingyu, 03/03/2020

subroutine calc_chemistry_new(iBlock)

  use ModSizeGitm
  use ModGITM
  use ModPlanet
  use ModRates
  use ModEUV
  use ModSources
  use ModInputs, only: iDebugLevel, UseIonChemistry, UseNeutralChemistry, &
                       UseNOChemiEquilibrium, UseN2AProfiles, f107, f107a
  use ModConstants
  use ModNOChemistry
  use ModTime, only: istep,utime,ijulianday
  use EUA_ModMsis90, only: meter6, gtd6
  use ModTime
  use ModIndicesInterfaces, only: get_HPI

  implicit none

  integer, intent(in) :: iBlock

  real :: IonSources(nIons), NeutralSources(nSpeciesTotal)
  real :: IonLosses(nIons), NeutralLosses(nSpeciesTotal)
  real :: DtSub, DtOld, DtTotal, DtMin, DtAve, Source, Reaction, tr, tr3, rr
  real :: te3, ti, tn, tn1, tn06, dtsubtmp, losstmp, dentmp, l, t,m1,m2,y1,y2,k1,k2
  real :: Ions(nIons), Neutrals(nSpeciesTotal)
  real :: tli(nIons), tsi(nIons), tln(nSpeciesTotal), tsn(nSpeciesTotal)

  integer :: iLon, iLat, iAlt, iIon, nIters, iNeutral, iVL, iError

  real :: lon, lat, alt
  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens
  real :: LonDeg, LatDeg, AltKm, LST
  real,dimension(7)    :: AP
  real :: hp

  real :: ChemicalHeatingSub,percent,o2ptotal
  real :: Emission(nEmissions), EmissionTotal(nEmissions)

  real, dimension(nLons,nLats,nAlts) :: &
       tr3d, tr33d,te12d, tr3m0443d, tr3m1163d, tr3m043d, &
       te3m0393d, te3m0693d, te33d, te3m053d,te3m073d,te12m0563d,te227d

  real, dimension(nLons,nLats,nAlts) :: &
       teffective_n2, teffective_o2, teffective_no, u2, mb, mbb, &
       k1_n2, k2_o2, k3_no

  real :: k1_n2_point, k2_o2_point, k3_no_point

  real :: te3m05,te3m07,te12m056, tr3m044, tr3m04, tr3m116
  real :: te3m039, te3m069, rr_opn2, te22m05
  real :: ionso, ionlo, neuso, neulo
  real :: NOP, NOL, NO_Emission


  logical :: UseNeutralConstituent(nSpeciesTotal)
  logical :: UseIonConstituent(nIons)
  !---------------------------------------------------------------------------

  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> start calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif


  UseNeutralConstituent = .true.
  UseIonConstituent     = .true.

  !  UseNeutralConstituent(iO_1D_) = .false.
  !  UseIonConstituent(iO_2PP_) = .false.
  !  UseIonConstituent(iO_2DP_) = .false.
  !
  !  UseNeutralConstituent(iN_4S_) = .false.
  !  UseNeutralConstituent(iN_2D_) = .false.
  !  UseNeutralConstituent(iO2_) = .false.
  !
  !  UseNeutralConstituent = .false.
  !  UseIonConstituent = .false.
  !  UseIonConstituent(1) = .true.
  !  UseIonConstituent(2) = .true.
  !  UseIonConstituent(3) = .true.
  !
  !  UseNeutralConstituent(1) = .true.
  !  UseNeutralConstituent(2) = .true.
  !  UseNeutralConstituent(3) = .true.

  !  open(unit=95,file='data.dat')
  DtMin = Dt

  if (.not.UseIonChemistry) return

  call report("Chemistry",2)
  call start_timing("calc_chemistry_new")

  DtAve = 0.0

  nIters=0

  !  AuroralIonRateS = 0.0

  u2 = IVelocity(1:nLons,1:nLats,1:nAlts,iEast_,iBlock)**2 + &
       IVelocity(1:nLons,1:nLats,1:nAlts,iNorth_,iBlock)**2 + &
       IVelocity(1:nLons,1:nLats,1:nAlts,iUp_,iBlock)**2

  mb  = 0.0
  mbb = 0.0

  ! This is from Schunk and Nagy 2nd Ed, formula 12.13 (pg 416)
  do iNeutral = 1, nSpeciesTotal
     ! Collisions should be better defined
     mbb = mbb + &
          (Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)) / &
          (mass(iNeutral) + MassI(iO_4SP_))
     mb  = mb + &
          (mass(iNeutral) * Collisions(1:nLons,1:nLats,1:nAlts,iVIN_)) / &
          (mass(iNeutral) + MassI(iO_4SP_))
  enddo

  mb = mb/mbb

  teffective_n2 = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) + Mass(iN2_)) * &
       (Mass(iN2_) + mb)/(3*Boltzmanns_Constant) * u2

  teffective_o2 = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) + Mass(iO2_)) * &
       (Mass(iO2_) + mb)/(3*Boltzmanns_Constant) * u2

  teffective_no = iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) + &
       MassI(iO_4SP_)/(MassI(iO_4SP_) + Mass(iNO_)) * &
       (Mass(iNO_) + mb)/(3*Boltzmanns_Constant) * u2

  where (teffective_n2 < 350.0)
     teffective_n2 = 350.0
  endwhere

  where (teffective_n2 > 6000.0)
     teffective_n2 = 6000.0
  endwhere

  where (teffective_o2 < 350.0)
     teffective_o2 = 350.0
  endwhere

  where (teffective_o2 > 350.0)
     teffective_o2 = 350.0
  endwhere

  where (teffective_no < 320.0)
     teffective_no = 320.0
  endwhere

  where (teffective_n2 <= 1700.0)
     k1_n2 =   1.533e-18 &
          - 5.920e-19*(teffective_n2/300.0) &
          + 8.600e-20*(teffective_n2/300.0)**2
  endwhere

  where (teffective_n2 > 1700.0)
     k1_n2 =   2.730e-18 &
          - 1.155e-18*(teffective_n2/300.0) &
          + 1.483e-19*(teffective_n2/300.0)**2
  endwhere

  k2_o2 =   2.820e-17 &
       - 7.740e-18*(teffective_o2/300.0) &
       + 1.073e-18*(teffective_o2/300.0)**2 &
       - 5.170e-20*(teffective_o2/300.0)**3 &
       + 9.650e-22*(teffective_o2/300.0)**4

  where (teffective_no <= 1500.0)
     k3_no =   8.360e-19 &
          - 2.020e-19*(teffective_no/300.0) &
          + 6.950e-20*(teffective_no/300.0)**2
  endwhere

  where (teffective_no > 1500.0)
     k3_no =   5.330e-19 &
          - 1.640e-20*(teffective_no/300.0) &
          + 4.720e-20*(teffective_no/300.0)**2 &
          - 7.050e-22*(teffective_no/300.0)**3
  endwhere

  tr3d = (iTemperature(1:nLons,1:nLats,1:nAlts,iBlock) &
       + Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
       TempUnit(1:nLons,1:nLats,1:nAlts)) / 2.0

  tr33d = tr3d/300.0
  te33d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/300.0
  te12d = eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)/1200.0
  te227d = -22740.0/eTemperature(1:nLons,1:nLats,1:nAlts,iBlock)

  te3m073d   = te33d**(-0.7)
  te12m0563d = te12d**(-0.56)
  te3m053d   = te33d**(-0.5)
  te3m0393d  = te33d**(-0.39)
  te3m0693d  = te33d**(-0.69) !J Yonker dissertation, Sheehan & StMaurice 2004 coef.

  tr3m0443d = tr33d**(-0.44)
  tr3m1163d  = tr33d**(-1.16)
  tr3m043d  = tr33d**(-0.4)

  m1 = ALOG(1.0/100000.0)/(115.0-100.0)
  k1 = 100000.0*exp(-m1*100.0)
  m2 = ALOG(1.0/1000.0)/(180.0-100.0)
  k2 = 1000.0*exp(-m2*100.0)


  do iLon = 1, nLons
     do iLat = 1, nLats
        do iAlt = 1, nAlts

           y1 = max(1.0,k1*exp(m1*altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0))
           y2 = max(1.0,k2*exp(m2*altitude_GB(iLon,iLat,iAlt,iBlock)/1000.0))
           NeutralSourcesTotal = 0.0
           NeutralLossesTotal = 0.0
           tr  = tr3d(iLon,iLat,iAlt)
           tr3 = tr33d(iLon,iLat,iAlt)
           te3 = te33d(iLon,iLat,iAlt)
           te3m05  = te3m053d(iLon,iLat,iAlt)
           te3m07  = te3m073d(iLon,iLat,iAlt)
           te3m069 = te3m0693d(iLon,iLat,iAlt)
           te3m039 = te3m0393d(iLon,iLat,iAlt)
           te12m056 = te12m0563d(iLon,iLat,iAlt)
           tr3m044 = tr3m0443d(iLon,iLat,iAlt)
           tr3m04  = tr3m043d(iLon,iLat,iAlt)
           te22m05  = eTemperature(iLon,iLat,iAlt,iBlock)**(.5) * &
                exp(te227d(iLon,iLat,iAlt))
           tr3m116  = tr3m1163d(iLon,iLat,iAlt)
           ti = iTemperature(iLon,iLat,iAlt,iBlock)
           tn = Temperature(iLon,iLat,iAlt,iBlock)*&
                TempUnit(iLon,iLat,iAlt)

           tn1 = exp(107.8/tn)
           tn06 = exp(67.5/tn)

           rr_opn2 = min(5.0e-19,4.5e-20*tr3**2)

           k1_n2_point = k1_n2(iLon,iLat,iAlt)
           k2_o2_point = k2_o2(iLon,iLat,iAlt)
           k3_no_point = k3_no(iLon,iLat,iAlt)

           DtTotal = 0.0
           EmissionTotal = 0.0

           Ions = IDensityS(iLon,iLat,iAlt,:,iBlock)

           Neutrals = NDensityS(iLon,iLat,iAlt,:,iBlock)

           niters = 0
           o2ptotal = 0

           ! NO stuff
           NOP = 0.0
           NOL = 0.0


           do while (DtTotal < Dt)

              ChemicalHeatingSub = 0.0
              ChemicalHeatingS = 0
              Emission = 0.0

              DtSub = Dt - DtTotal

              IonSources = 0.0
              NeutralSources = 0.0
              IonLosses  = 0.0
              NeutralLosses = 0.0

              ! ----------------------------------------------------------
              ! O2 -> 2O
              ! ----------------------------------------------------------
              rr=EuvDissRateS(iLon,iLat,iAlt,iO2_,iBlock)

              Reaction = rr * &
                   Neutrals(iO2_)

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 2*Reaction

              ! -----------------------------------------------------------
              ! O+O+M -> O2+M
              ! -----------------------------------------------------------
              rr=9.59e-34 * exp(480./tn)
              rr= rr*1.e-12  !cm6s-1-->m6s-1

              Reaction = rr * Neutrals(iO_3P_)**2 *&
                   (Neutrals(iO2_)+ &
                   Neutrals(iO_3P_)+ &
                   Neutrals(iN2_))

              NeutralLosses(iO_3P_) = NeutralLosses(iO_3P_) + 2*Reaction
              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction

              ! -----------------------------------------------------------
              ! O+O2+M -> O3+M
              ! -----------------------------------------------------------
!!$              rr=6.0e-34 * exp(300./tn)
!!$              rr=rr*1.e-12  !cm6s-1-->m6s-1
!!$
!!$              Reaction = rr * Neutrals(iO_3P_) *&
!!$                    Neutrals(iO2_) * &
!!$                   (Neutrals(iO2_) + &
!!$                    Neutrals(iO_3P_)  + &!$                    Neutrals(iN2_))
!!$
!!$              NeutralLosses(iO_3P_) = NeutralLosses(iO_3P_) + Reaction
!!$              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction


              ! ----------------------------------------------------------
              ! N2(A)
              ! - added by C Lin 2016
              ! ----------------------------------------------------------                                    

              ! -----------------                                                                             
              ! What's temporarily missing here is the N2A production rate                                    
              ! -----------------                                                                             
              ! needs to be added here                                                                        
              !                                                                                               


              !
              ! Loss of N2(A) - added by C Lin 2016                                                           
              !
              ! N2(A) + O(3P) -> NO + N(2D)
              ! rr = 1e-19 * sqrt(T/298) for NOx1d SNOE
              !      [Dilecce and DeBenedictis, 1999]
              ! no T dep. for NOx1D EUVAC [Hill et al., 2000]
              !                                                                                               

              if (UseN2AProfiles) then
                  rr = 1e-19 * sqrt(tn/298.0)
                  Reaction = rr * &
                             Neutrals(iN2_A_) * &
                             Neutrals(iO_3P_)

                  NeutralLosses(iN2_A_)  = NeutralLosses(iN2_A_) + Reaction
                  NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_) + Reaction
                  NeutralSources(iNO_)   = NeutralSources(iNO_) + Reaction
                  NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction


                  !if (mod(utime, 10.0) == 0.0) then
                  !    write(*,*) 'nO_3P = ', NeutralLosses(iO_3P_)
                  !    write(*,*) 'nN2_A = ', NeutralLosses(iN2_A_)
                  !    write(*,*) 'nNO = ', NeutralLosses(iNO_)
                  !    write(*,*) 'nN_2D = ', NeutralLosses(iN_2D_)
                  !endif

              endif



              ! ----------------------------------------------------------
              ! N2 -> 2N
              ! ----------------------------------------------------------

              Reaction = EuvDissRateS(iLon,iLat,iAlt,iN2_,iBlock) * &
                   Neutrals(iN2_)

              !              Reaction = 5.0e-12 * &
              !                   Neutrals(iN2_)

              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + .25*Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + .60*Reaction

              ! Solar EUV

              ! ----------------------------------------------------------
              ! N2+
              ! ----------------------------------------------------------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock) * &
                   Neutrals(iN2_)
!              if (UseTCIParameterizationS) then
!                 if (ScalingEUVIonRate .eq. 1.0) then
!                      call get_HPI(CurrentTime, HP, iError)
!                      ScalingEUVIonRate = 1.0 / (1.50137e-2 * (hp - 4.67) + 2.91889)
!                  endif
!                  write(*,*) "ScalingEUVIonRate in Chem = ", ScalingEUVIonRate
!                  Reaction = ScalingEuvIonRate * EuvIonRateS(iLon,iLat,iAlt,iN2P_,iBlock) * &
!                             Neutrals(iN2_)
!              endif

              IonSources(iN2P_)   = IonSources(iN2P_)   + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! Aurora

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock) + &
                   IonPrecipIonRateS(iLon,iLat,iAlt,iN2_, iBlock)
              if (UseTCIParameterizationGM) then
                  ScalingIonPrecipIonRate = 1.0 / (8.48315e-3 * (f107 - 70.0) + 0.915558)
                  Reaction = AuroralIonRateS(iLon,iLat,iAlt,iN2_, iBlock) + &
                              ScalingIonPrecipIonRate * IonPrecipIonRateS(iLon,iLat,iAlt,iN2_, iBlock)
              endif

              IonSources(iN2P_)   = IonSources(iN2P_) + Reaction
              NeutralLosses(iN2_) = NeutralLosses(iN2_) + Reaction

              ! O+(2D) + N2 -> N2+ + O + 1.33 eV

              Reaction = &
                   8.0e-16 * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN2_)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
              IonLosses(iO_2DP_)   = IonLosses(iO_2DP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 1.33

              ChemicalHeatingS(iop2d_n2) =  &
                   ChemicalHeatingS(iop2d_n2) + &
                   Reaction * 1.33

              ! O+(2P) + N2 -> N2+ + O + 3.02 eV

              Reaction = &
                   4.8e-16 * &
                   Ions(iO_2PP_) * &
                   Neutrals(iN2_)

              IonSources(iN2P_)    = IonSources(iN2P_)   + Reaction
              NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + Reaction
              IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)  + Reaction
              NeutralLosses(iN2_)  = NeutralLosses(iN2_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 3.02

              ChemicalHeatingS(iop2p_n2) =  &
                   ChemicalHeatingS(iop2p_n2) + &
                   Reaction * 3.02

              ! N2+ + O2 -> O2+ + N2 + 3.53 eV
              ! rr = 5.10e-17 * t  [Fox & Sung 2001]  - C Lin 2015
              ! so 2016 code has ti3m116, I used tr3m116, which one is right?
              ! I'm keeping it ti3m116 for now

              rr = 5.10e-17 * tr3m116

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iO2_)

              IonSources(iO2P_)    = IonSources(iO2P_)   + Reaction
              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)  + Reaction
              NeutralLosses(iO2_)  = NeutralLosses(iO2_) + Reaction

              o2ptotal = o2ptotal + reaction
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + Reaction * 3.53

              ChemicalHeatingS(in2p_o2) =  &
                   ChemicalHeatingS(in2p_o2) + &
                   Reaction * 3.53

              ! N2+ + O -> NO+ + N(2D) + 0.70 eV
              !!!!!          -> NO+ + N(4S) + 3.08 eV
              !
              ! N2+ + O(3P) --0.90--> N(2D) + NO+
              !         --0.05--> N(4S) + NO+
              !         --0.05--> O+(4S) + N2 + 1.96eV
              ! Section 3.1.1 of Yonker 2013; Scott et al. 1999
              ! rr = 1.4e-16 * tr3m044  Fox & Sung 2001
              ! -- modified by Cissi Lin 2015 --
              ! again, keep it the 2012 rate variable (ti3m044 in v2016)
              ! -- C Lin 2016
              rr = 1.4e-16 * tr3m044

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iO_3P_)

              IonSources(iNOP_)      = IonSources(iNOP_)      + 0.95*Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 0.90*Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 0.05*Reaction
              IONSources(iO_4SP_) = IonSources(iO_4SP_) + 0.05*Reaction
              NeutralSources(iN2_) = NeutralSources(iN2_) + 0.05*Reaction
              IonLosses(iN2P_)       = IonLosses(iN2P_)       + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   0.90 * Reaction * 0.70

              ChemicalHeatingS(in2p_o) =  &
                   ChemicalHeatingS(in2p_o) + &
                   0.90 * Reaction * 0.70

              ! Heating from:
              ! N2+ + O -> O+(4S) + N2 + 1.96 eV
              ! Cissi Lin
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   0.05 * Reaction * 1.96

              ChemicalHeatingS(in2p_o) =  &
                   ChemicalHeatingS(in2p_o) + &
                   0.05 * Reaction * 1.96



              ! N2+ + e -> 2 N(2D) + 1.04 eV
              !
              ! -------
              ! N2+ + e --0.52--> N(2D) + N(2D)
              !         --0.37--> N(4S) + N(2D)
              !         --0.11--> N(4S) + N(2P)  [Peterson et al 1998]
              ! rr = 2.2e-13 * te3m039  [Sheehan & StMaurice 2004]
              ! -- modified by Cissi Lin 2015 --

              rr = 2.2e-13 * te3m039

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Ions(ie_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 1.41*Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 0.48*Reaction
              NeutralSources(iN_2P_) = NeutralSources(iN_2P_) + 0.11*Reaction
              IonLosses(iN2P_)       = IonLosses(iN2P_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   0.52 * Reaction * 1.04

              ChemicalHeatingS(in2p_e) =  &
                   ChemicalHeatingS(in2p_e) + &
                   0.52 * Reaction * 1.04


              ! N2+ + NO -> NO+ + N2 + 6.33 eV

              !rr = 3.3e-16
              rr = 3.6e-16  ! Richards - changed from 2012 value (3.3e-16)

              Reaction = &
                   rr * &
                   Ions(iN2P_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_) = NeutralSources(iN2_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)     + Reaction
              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonLosses(iN2P_)     = IonLosses(iN2P_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 6.33


              ! ----------------------------------------------------------
              ! O2+
              ! ----------------------------------------------------------

              ! -----------
              ! Solar EUV
              ! -----------

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO2P_,iBlock) * &
                   Neutrals(iO2_)

              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              o2ptotal = o2ptotal + reaction

              ! -----------
              ! Aurora
              ! -----------

              Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO2_, iBlock) + &
                   IonPrecipIonRateS(iLon,iLat,iAlt,iO2_, iBlock)

              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction

              o2ptotal = o2ptotal + reaction

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction

              ! -----------
              ! O+(4S) + O2 -> O2+ + O + 1.55 eV
              ! -----------

              !! rr = 2.0e-17 * tr3m04

              rr = k2_o2_point

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              o2ptotal = o2ptotal + reaction

              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction
              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.55

              ChemicalHeatingS(iop_o2) =  &
                   ChemicalHeatingS(iop_o2) + &
                   Reaction * 1.55

              ! -----------
              ! O+(2D) + O2 -> O2+ + O + 4.865 eV
              ! -----------

              rr = 7.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iO2_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iO2P_)   = IonSources(iO2P_)   + Reaction
              o2ptotal = o2ptotal + reaction
              NeutralLosses(iO2_) = NeutralLosses(iO2_) + Reaction/1.0
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_) + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.865

              ! -----------
              ! N+ + O2 -> O2+ + N(4S) + 2.486 eV
              ! -----------

              rr = 1.1e-16

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iO2P_)      = IonSources(iO2P_)      + Reaction
              o2ptotal = o2ptotal + reaction
              NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)       + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.486

              ChemicalHeatingS(inp_o2) =  &
                   ChemicalHeatingS(inp_o2) + &
                   Reaction * 2.486

              ! -----------                                                                       
              ! N+ + O2 -> O2+ + N(2D) + 0.1 eV                                                   
              ! -----------                                                                       
              !                                                                                   
              ! C Lin: BR issue here, products are everywhere                                     
              ! N+ + O2 -0.50-> O2+ + N(2D) + 0.1 eV                                              
              !         -0.42-> NO+ + O(3P) + 6.67 eV                                             
              !                 NO+ + O(1D) + 4.71 eV                                             
              !         -0.08-> O+(4S) + NO + 2.31 eV                                             
              ! rr = 5.5e-16  [Midey et al. 2006]                                                 

              !if (ti<=1000.0) then                                                               
              !    rr =  0.825e-16 * ti3m045                                                      
              !else                                                                               
              !    rr = 1.425e-16                                                                 
              !endif                                                                              

              !rr = 2.0e-16
              rr = 5.5e-16

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 0.50 * Reaction
              IonSources(iO2P_)      = IonSources(iO2P_)      + 0.50 * Reaction
              o2ptotal = o2ptotal + reaction

              IonSources(iNOP_)      = IonSources(iNOP_)      + 0.42 * Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.42 * Reaction

              IonSources(iO_4SP_)    = IonSources(iO_4SP_)    + 0.08 * Reaction
              NeutralSources(iNO_)   = NeutralSources(iNO_)   + 0.08 * Reaction

              NeutralLosses(iO2_)    = NeutralLosses(iO2_)    + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction


              ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  0.50 * Reaction * 0.07 + &
                  0.42 * Reaction * 4.44 + &
                  0.08 * Reaction * 1.54

              ChemicalHeatingS(inp_o2) =  &
                  ChemicalHeatingS(inp_o2) + &
                  0.50 * Reaction * 0.1


              ! -----------
              ! O2+ + e -> O(1D) + O(1D) + 3.06 eV
              ! O2+ + e -> O(3P) + O(1D) + 5.02 eV
              ! O2+ + e -> O(3P) + O(3P) + 6.99 eV
              ! -----------

              !              rr = 1.9e-13 * te3m05
              !rr = 2.4e-13 * te3m07  ! commented out by C Lin
              rr = 1.95e-13 * te3m07  ![Sheehan and St-Maurice, 2004]
              !              if (eTemperature(iLon,iLat,iAlt,iBlock) .lt. 1200)
              !                 rr = 10.9e-13* te12m056
              !              else
              !                 rr = 7.4e-14* te12m056
              !              endif

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Ions(ie_)

              if (UseNeutralConstituent(iO_1D_)) then
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.22*Reaction * 2.0
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.42*Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.42*Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.31*Reaction * 2.0
                 ! This really should be 0.05 to O(1D) and O(1S)
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.05*Reaction
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.05*Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 6.99 * 0.22 + &
                      Reaction * 5.02 * 0.42 + &
                      Reaction * 3.06 * 0.31

                 ChemicalHeatingS(io2p_e) = &
                      ChemicalHeatingS(io2p_e) + &
                      Reaction * 6.99 * 0.22 + &
                      Reaction * 5.02 * 0.42 + &
                      Reaction * 3.06 * 0.31

              else

                 NeutralSources(iO_3P_)    = NeutralSources(iO_3P_) + Reaction * 2.0

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 5.0

                 ChemicalHeatingS(io2p_e) = &
                      ChemicalHeatingS(io2p_e) + &
                      Reaction * 5.0

              endif

              IonLosses(iO2P_)      = IonLosses(iO2P_)   + Reaction

              ! -----------
              ! O2+ + N(4S) -> NO+ + O + 4.21 eV
              ! -----------
              ! energy release value chaged from 2012 verison (4.25 eV)
              !rr = 1.8e-16
              rr = 1.0e-16 ! Richards 2011 [Scott et al., 1999]                                   


              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNOP_)     = IonSources(iNOP_)     + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.

              ChemicalHeatingS(io2p_n) =  &
                   ChemicalHeatingS(io2p_n) + &
                   Reaction * 4.25

              ! -----------
              ! O2+ + NO -> NO+ + O2 + 2.813 eV
              ! -----------

              !rr = 4.4e-16 for v2012
              rr = 4.5e-16 ! schunk and nagy                                                      

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iNO_)

              NeutralSources(iO2_) = NeutralSources(iO2_) + Reaction
              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
              NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.813

              ChemicalHeatingS(io2p_no) =  &
                   ChemicalHeatingS(io2p_no) + &
                   Reaction * 2.813

!--- Comment'd out C Lin
!
!              ! -----------
!              ! O2+ + N2 -> NO+ + NO + 0.9333 eV
!              ! -----------
!
!              rr = 5.0e-22
!
!              Reaction = &
!                   rr * &
!                   Ions(iO2P_) * &
!                   Neutrals(iN2_)
!
!              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!
!
!              IonSources(iNOP_)    = IonSources(iNOP_)    + Reaction
!              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
!              IonLosses(iO2P_)     = IonLosses(iO2P_)     + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.9333
!
!              ChemicalHeatingS(io2p_n2) =  &
!                   ChemicalHeatingS(io2p_n2) + &
!                   Reaction * 0.9333
!

              ! ----------------------------------------------------------
              ! O(4S)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_4SP_,iBlock) * &
                   Neutrals(iO_3P_)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              ! Aurora goes 0.4, 0.4, 0.2 into O(4S), O(2D) and O(2P) respectively
              Reaction = 0.4 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.4 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora goes 0.4, 0.4, 0.2 into O(3P), O(2D) and O(2P) respectively
              Reaction = 0.4 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.4 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora goes 0.4, 0.4, 0.2 into O(3P), O(2D) and O(2P) respectively
              Reaction = 0.2 * AuroralIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock) + &
                   0.2 * IonPrecipIonRateS(iLon,iLat,iAlt,iO_3P_, iBlock)

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2D) + O -> O+(4S) + O(3P) + 3.31 eV
              ! O+(2D) + O -> O+(4S) + O(1D) + 1.35 eV
              ! -----------

              rr = 1.0e-17

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iO_3P_)

              ! We create and loose the same amount of O (when O(1D) is not used...

              if (UseNeutralConstituent(iO_1D_)) then

                 NeutralSources(iO_3P_)  = NeutralSources(iO_3P_) + 0.5 * Reaction
                 NeutralSources(iO_1D_)  = NeutralSources(iO_1D_) + 0.5 * Reaction

                 NeutralLosses(iO_3P_)   = NeutralSources(iO_3P_)  + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.31 * 0.5 + &
                      Reaction * 1.35 * 0.5
              else

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.31

              endif

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ! -----------
              ! O+(2D) + e -> O+(4S) + e + 3.31 eV
              ! -----------

              rr = 7.8e-14 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.31

              ChemicalHeatingS(iop2d_e) =  &
                   ChemicalHeatingS(iop2d_e) + &
                   Reaction * 3.31

!              ! -----------
!              ! O+(2D) + N2 -> O+(4S) + N2 + 3.31 eV
!              ! -----------
!
!              rr = 8.0e-16
!
!              Reaction = &
!                   rr * &
!                   Ions(iO_2DP_) * &
!                   Neutrals(iN2_)
!
!              ! We create and loose the same amount of N2
!              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
!              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 3.31
!
!              ChemicalHeatingS(iop2d_n2) =  &
!                   ChemicalHeatingS(iop2d_n2) + &
!                   Reaction * 3.31


              ! -----------
              ! O+(2P) + O -> O+(4S) + O + 5.0 eV
              ! -----------

              !rr = 5.2e-17 (v2012)
              rr = 4.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iO_3P_)

              ! We create and loose the same amount of O
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.0

              ChemicalHeatingS(iop2p_o) =  &
                   ChemicalHeatingS(iop2p_o) + &
                   Reaction * 5.0

              ! -----------
              ! O+(2P) + e -> O+(4S) + e + 5.0 eV
              ! -----------

              rr = 4.0e-14 * te3m05

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.0

              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction * 5.0

              ! -----------
              ! O+(2P) -> O+(4S) + 2470A
              ! -----------

              rr = 0.047

              Reaction = &
                   rr * &
                   Ions(iO_2PP_)

              ! We create and loose the same amount of e
              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE2470_) = Emission(iE2470_) + Reaction

              !              ! -----------
              !              ! H+ + O -> O+(4S) + H + 0.0 eV
              !              ! -----------
              !
              !              rr = 6.0e-10/1.0e6 * (8.0/9.0) * sqrt( &
              !                   (ti + tn/4) * (tn + ti/16))
              !
              !              Reaction = &
              !                   rr * &
              !                   Ions(iHP_) * &
              !                   Neutrals(iO_3P_)
              !
              !              NeutralSources(iH_) = NeutralSources(iH_) + Reaction
              !              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              !              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction
              !              IonLosses(iHP_)     = IonLosses(iHP_)     + Reaction
              !
              !              ChemicalHeatingSub = &
              !                   ChemicalHeatingSub + &
              !                   Reaction * 0.0


              ! -----------
              ! O+(4S) + N2 -> NO+ + N(4S) + 1.10 eV
              ! -----------

              !!              rr = rr_opn2
              rr = k1_n2_point

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iN2_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
              IonLosses(iO_4SP_)     = IonLosses(iO_4SP_)     + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.10

              ChemicalHeatingS(iop_n2) =  &
                   ChemicalHeatingS(iop_n2) + &
                   Reaction * 1.10

              ! -----------
              ! O+(4S) + NO -> NO+ + O + 4.36 eV
              ! -----------

              !!              rr = 8.0e-19
              rr = k3_no_point

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iNO_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction
              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 4.36

              !              ! -----------
              !              ! O+(4S) + H -> H+ + O + 0.0 eV
              !              ! -----------
              !
              !              rr = 6.0e-10/1.0e6
              !
              !              Reaction = &
              !                   rr * &
              !                   Ions(iO_4SP_) * &
              !                   Neutrals(iH_)
              !
              !              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              !              IonSources(iHP_)    = IonSources(iHP_)    + Reaction
              !              NeutralLosses(iH_)  = NeutralLosses(iH_)  + Reaction
              !              IonLosses(iO_4SP_)  = IonLosses(iO_4SP_)  + Reaction
              !
              !              ChemicalHeatingSub = &
              !                   ChemicalHeatingSub + &
              !                   Reaction * 0.0

              ! -----------
              ! O+(4S) + N(2D) -> N+ + O + 1.45 eV
              ! -----------

              rr = 1.3e-16

              Reaction = &
                   rr * &
                   Ions(iO_4SP_) * &
                   Neutrals(iN_2D_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              IonLosses(iO_4SP_)    = IonLosses(iO_4SP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.45

              ! ----------------------------------------------------------
              ! O(2D)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2DP_,iBlock) * &
                   Neutrals(iO_3P_)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              !Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_, iBlock)

              !IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              !NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! -----------
              ! O+(2P) + e -> O+(2D) + e + 1.69 eV
              ! -----------

              !rr = 1.5e-13 * te3m05 (v2012)
              rr = 1.84e-13 * te3m05 !(v2016)

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.69

              ChemicalHeatingS(iop2p_e) =  &
                   ChemicalHeatingS(iop2p_e) + &
                   Reaction * 1.69


              ! -----------
              ! O+(2P) -> O+(2D) + 7320A
              ! -----------

              rr = 0.171

              Reaction = &
                   rr * &
                   Ions(iO_2PP_)

              IonSources(iO_2DP_) = IonSources(iO_2DP_) + Reaction
              IonLosses(iO_2PP_)  = IonLosses(iO_2PP_)  + Reaction

              Emission(iE7320_) = Emission(iE7320_) + Reaction


              ! -----------
              ! O+(2D) + N2 -> NO+ + N + 4.41 eV
              ! -----------                                                                       

              rr = 2.5e-17

                    Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN2_)

              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
              IonSources(iNOP_) = IonSources(iNOP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
              NeutralSources(iN_4S_)  = NeutralSources(iN_4S_)  + Reaction


              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.007


              ! -----------
              ! O+(2D) + NO -> NO+ + O + 4.37 eV
              ! -----------

             rr = 1.2e-15

             Reaction = &
                  rr * &
                  Ions(iO_2DP_) * &
                  Neutrals(iNO_)

             NeutralLosses(iNO_)  = NeutralLosses(iNO_)  + Reaction
             IonSources(iNOP_) = IonSources(iNOP_) + Reaction
             IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction
             NeutralSources(iO_3P_)  = NeutralSources(iO_3P_)  + Reaction


             ChemicalHeatingSub = &
                  ChemicalHeatingSub + &
                  Reaction * 2.85


              ! -----------
              ! O+(2D) -> O+(4S) + 3726A
              ! -----------

              rr = 7.7e-5

              Reaction = &
                   rr * &
                   Ions(iO_2DP_)

              IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              IonLosses(iO_2DP_)  = IonLosses(iO_2DP_)  + Reaction

              Emission(iE3726_) = Emission(iE3726_) + Reaction

              ! ----------------------------------------------------------
              ! O(2P)+
              ! ----------------------------------------------------------

              ! Solar EUV

              Reaction = EuvIonRateS(iLon,iLat,iAlt,iO_2PP_,iBlock) * &
                   Neutrals(iO_3P_)

              IonSources(iO_2PP_) = IonSources(iO_2PP_) + Reaction
              NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction

              ! Aurora

              !Reaction = AuroralIonRateS(iLon,iLat,iAlt,iO_, iBlock)

              !IonSources(iO_4SP_) = IonSources(iO_4SP_) + Reaction
              !NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_)  + Reaction


! --- Comment'd out by C Lin
!
!              ! -----------
!              ! O+(2P) + N2 -> N+ + NO + 0.70 eV
!              ! -----------
!
!              rr = 1.0e-16
!
!              Reaction = &
!                   rr * &
!                   Ions(iO_2PP_) * &
!                   Neutrals(iN2_)
!
!              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction
!              IonSources(iNP_)     = IonSources(iNP_)     + Reaction
!              NeutralLosses(iN2_)  = NeutralLosses(iN2_)  + Reaction
!              IonLosses(iO_2PP_)   = IonLosses(iO_2PP_)   + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 0.70
!
!              ChemicalHeatingS(iop2p_n2) =  &
!                   ChemicalHeatingS(iop2p_n2) + &
!                   Reaction * 0.70


              ! ----------------------------------------------------------
              ! N+
              ! ----------------------------------------------------------

              ! -----------
              ! O2+ + N(2D) -> N+ + O2 + 0.0 eV
              ! -----------
              ! Cissi Lin: I don't know why this reaction was comment'd out in 2016 version
              !        No reason provided. So I resumed it with Richards reactions & values.
              !rr = 8.65e-17  ! C Lin: I don't where this rate is coming from
              !rr = 2.5e-16  ! Rees (v2012)
              rr = 8.65e-17  ! Richards 2011

              rr = 2.5e-16

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iN_2D_)

              NeutralSources(iO2_)  = NeutralSources(iO2_)  + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              IonLosses(iO2P_)      = IonLosses(iO2P_)      + Reaction
              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.0



              ! Cissi Lin: newly added                                                            
              ! -----------                                                                       
              ! O2+ + N(2D) -> NO+ + O                                                            
              ! -----------                                                                       

              rr = 1.8e-16  ! Richards 2011                                                       

              Reaction = &
                 rr * &
                 Ions(iO2P_) * &
                 Neutrals(iN_2D_)

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonSources(iNOP_)      = IonSources(iNOP_)      + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction
              IonLosses(iO2P_)       = IonLosses(iO2P_)       + Reaction


              !              ! -----------
              !              ! He+ + N2 -> N+ + N + He + 0.28 eV
              !              ! -----------
              !
              !              rr = 1.2e-9/1.0e6
              !
              !              Reaction = &
              !                   rr * &
              !                   Ions(iHeP_) * &
              !                   Neutrals(iN2_)
              !
              !              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              !              NeutralSources(iHe_)   = NeutralSources(iHe_)   + Reaction
              !              IonSources(iNP_)       = IonSources(iNP_)       + Reaction
              !              NeutralLosses(iN2_)    = NeutralLosses(iN2_)    + Reaction
              !              IonLosses(iHeP_)       = IonLosses(iHeP_)       + Reaction
              !
              !              ChemicalHeatingSub = &
              !                   ChemicalHeatingSub + &
              !                   Reaction * 0.28

!              ! -----------
!              ! O+(2P) + N -> N+ + O + 2.7 eV
!              ! -----------
!
!              rr = 1.0e-16
!
!              Reaction = &
!                   rr * &
!                   Ions(iO_2PP_) * &
!                   Neutrals(iN_4S_)
!
!              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
!              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
!              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
!              IonLosses(iO_2PP_)    = IonLosses(iO_2PP_)    + Reaction
!
!              ChemicalHeatingSub = &
!                   ChemicalHeatingSub + &
!                   Reaction * 2.7


              ! -----------
              ! N+ + NO --> NO+ + N(4S) + 3.4 eV
              ! -----------
              !
              ! N+ + NO -0.91-> NO+ + N(4S) + 3.4 eV
              !         -0.07-> N2+ + O + 2.2 eV
              !         -0.02-> N2 + O+(4S)

              !rr = 4.72e-16 * ti3m024
              ! Cissi Lin: reaction rate from Midey 2004
              rr = 6.5e-15 * ti ** (-0.44)

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iNO_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 0.91 * Reaction
              IonSources(iNOP_)      = IonSources(iNOP_)      + 0.91 * Reaction

              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.07 * Reaction
              IonSources(iN2P_)      = IonSources(iN2P_)      + 0.07 * Reaction

              NeutralSources(iN2_)   = NeutralSources(iN2_)   + 0.02 * Reaction
              IonSources(iO_4SP_)    = IonSources(iO_4SP_)    + 0.02 * Reaction

              NeutralLosses(iNO_)    = NeutralLosses(iNO_)    + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   0.91 * Reaction * 2.267 + &
                   0.07* Reaction * 1.47



              ! -----------
              ! O+(2D) + N -> N+ + O + 1.0 eV
              ! -----------

              !rr = 7.5e-17 (v2012)
              rr = 1.5e-16 !(v2016)

              Reaction = &
                   rr * &
                   Ions(iO_2DP_) * &
                   Neutrals(iN_4S_)

              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              IonSources(iNP_)      = IonSources(iNP_)      + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              IonLosses(iO_2DP_)    = IonLosses(iO_2DP_)    + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.0


! Cissi Lin: This reaction wasn't written clearly.                                                
              ! -----------
              ! N+ + O(3P) -> O+ + N(4S) + 0.93 eV
              ! -----------

              !rr = 2.2e-18
              ! C Lin: changed to the ACE1D rate (4.5e-18)
              !        Table D4 in Yonker 2013
              rr = 4.5e-18

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO_3P_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              IonSources(iO_4SP_)    = IonSources(iO_4SP_)    + Reaction
              NeutralLosses(iO_3P_)     = NeutralLosses(iO_3P_)     + Reaction
              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 0.496

              ChemicalHeatingS(inp_o) =  &
                   ChemicalHeatingS(inp_o) + &
                   Reaction * 0.93


              !              ! -----------
              !              ! N+ + H -> H+ + N + 0.90 eV
              !              ! -----------
              !
              !              rr = 3.6e-12/1.0e6
              !
              !              Reaction = &
              !                   rr * &
              !                   Ions(iNP_) * &
              !                   Neutrals(iH_)
              !
              !              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              !              IonSources(iHP_)       = IonSources(iHP_)       + Reaction
              !              NeutralLosses(iH_)     = NeutralLosses(iH_)     + Reaction
              !              IonLosses(iNP_)        = IonLosses(iNP_)        + Reaction
              !
              !              ChemicalHeatingSub = &
              !                   ChemicalHeatingSub + &
              !                   Reaction * 0.9

              ! ----------------------------------------------------------
              ! NO+
              ! ----------------------------------------------------------

              ! -----------
              !!!!!!! NO+ + e -> O + N(4S) + 2.75 eV  (0.22)
              ! NO+ + e -> O + N(2D) + 0.38 eV  (0.78)
              !
              ! NO+ + e --0.95--> N(2D) + O(3P)  (0.95)
              ! NO+ + e --0.05--> N(4S) + O  (0.05)  Hellberg et al 2003
              ! rr = 3.5e-13 * te3m069  Sheehan & StMaurice 2004
              ! Cissi Lin 2015
              ! -----------

              rr = 3.5e-13 * te3m069  !modified from earlier 4.2e-13 * te3m085

              Reaction = &
                   rr * &
                   Ions(iNOP_) * &
                   Ions(ie_)


              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 0.95*Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 0.05*Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
              IonLosses(iNOP_)       = IonLosses(iNOP_)       +      Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   0.95 * Reaction * 0.38 + &
                   0.05 * Reaction * 2.77

              ChemicalHeatingS(inop_e) =  &
                   ChemicalHeatingS(inop_e) + &
                   0.95 * Reaction * 0.38 + &
                   0.05 * Reaction * 2.77


              ! ----------------------------------------------------------
              ! N(4S)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + e -> N(4S) + e + 2.38 eV
              ! -----------

              rr = 5.5e-16 * te3 ** (0.5)

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Ions(ie_)

              ! We create and loose the same amount of e
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 2.38

              ChemicalHeatingS(in2d_e) =  &
                   ChemicalHeatingS(in2d_e) + &
                   Reaction * 2.38

              ! -----------
              ! N(2D) + O -> N(4S) + O(3P) + 2.38 eV
              ! N(2D) + O -> N(4S) + O(1D) + 0.42 eV
              ! -----------
              !rr = 2.0e-18 (v2012)
              ! Cissi Lin: rate from Table D6 of Yonker 2013                                      
              rr = 1.65e-19 * exp(-260.0/tn)

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO_3P_)

              if (UseNeutralConstituent(iO_1D_)) then
                 NeutralSources(iO_3P_)    = NeutralSources(iO_3P_)     + 0.9 * Reaction
                 NeutralSources(iO_1D_)    = NeutralSources(iO_1D_)  + 0.1 * Reaction
                 NeutralLosses(iO_3P_)     = NeutralLosses(iO_3P_)  + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      0.9 * Reaction * 2.38 + &
                      0.1 * Reaction * 0.42

                 ChemicalHeatingS(in2d_o) =  &
                      ChemicalHeatingS(in2d_o) + &
                      0.9 * Reaction * 2.38 + &
                      0.1 * Reaction * 0.42

              else

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 2.38

              endif

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction


              ! Cissi Lin: newly added
              ! -----------
              ! N(2D) + N2 -> N(4S) + N2
              ! -----------
              rr = 1.74e-20
              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                              Neutrals(iN2_)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iN_4S_)   = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iN2_)   = NeutralLosses(iN2_)   + Reaction



              ! -----------
              ! N(2D) -> N(4S) + 5200A
              ! -----------

              rr = 1.06e-5

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2D_)  = NeutralLosses(iN_2D_)  + Reaction
              Emission(iE5200_) = Emission(iE5200_) + Reaction

! Cissi Lin: This reaction is not in ACE1D, comment it out
!              ! -----------
!              ! NO -> N(4S) + O
!              ! -----------
!
!              !              rr = 8.3e-6
!              rr=4.5e-6*exp(-1.e-8*(Neutrals(iO2_)*1.e-6)**0.38)
!
!              Reaction = &
!                   rr * &
!                   Neutrals(iNO_)
!
!              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
!              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
!              NeutralLosses(iNO_)    = NeutralLosses(iNO_)    + Reaction


              
              ! -----------
              ! N(4S) + O2 -> NO + O + 1.385 eV
              ! -----------

              !rr = 4.4e-18 * exp(-3220/tn)  ! commented out by C Lin
              ! Cissi Lin: C16 in Sander 2015
              rr = 1.5e-17 * exp(-3600/tn)  ! [Sultanou et al., 2006]

              Reaction = &
                   rr * &
                   Neutrals(iN_4S_) * &
                   Neutrals(iO2_)

              NeutralSources(iNO_)  = NeutralSources(iNO_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 1.385

              ChemicalHeatingS(in_o2) =  &
                   ChemicalHeatingS(in_o2) + &
                   Reaction * 1.385



              ! Emission caused by Higher Vibrational State from N(4S) + O2 Reaction
              if (UseNOHigherVibrationalStates) then

                  do iVL = 1, 10

                      Emission(i5300_) = Emission(i5300_) + &
                         iVL * Planck_Constant * Speed_Light / 5.3e-6 * Adv1(iVL) / &
                         (Adv1(iVL) + Adv2(iVL) + VSP_f * VSP_k(iVL) * Neutrals(iO_3P_)) &
                         * VSP_g4s(iVL+1) * Reaction

                  enddo

              endif





              ! -----------
              ! N(4S) + NO -> N2 + O + 3.25 eV
              ! -----------
              !rr = 1.5e-18 * sqrt(tn)  ! commented out by C Lin
              ! rr = 3.4e-17
              rr = 2.1e-17 * exp(100./tn)

              Reaction = &
                   rr * &
                   Neutrals(iN_4S_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_4S_) = NeutralLosses(iN_4S_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 3.25

              ChemicalHeatingS(ino_n) =  &
                   ChemicalHeatingS(ino_n) + &
                   Reaction * 3.25


              ! -----------
              ! N(2P) -> N(2D) + 10400A
              ! -----------

              rr = 7.9e-2

              Reaction = &
                   rr * &
                   Neutrals(iN_2P_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction

              Emission(iE10400_) = Emission(iE10400_) + Reaction

              ! Cissi Lin: newly added                                                            
              ! -----------                                                                       
              ! N(2P) + e -> N(2D) + e                                                            
              ! -----------                                                                       
              rr = 9.5e-15
              Reaction = &
                    rr * &
                    Neutrals(iN_2P_) * &
                    Ions(ie_)

              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction
              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction


              ! Cissi Lin: newly added                                                            
              ! -----------                                                                       
              ! N(2P) + e -> N(4S) + e                                                            
              ! -----------                                                                       
              rr = 1.6e-18 * te3 ** (0.85)
              Reaction = &
                    rr * &
                    Neutrals(iN_2P_) * &
                                Ions(ie_)

              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + Reaction
              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction


              ! -----------                                                                       
              ! N(2P) + O(3P) -0.50-> NO+ + e                                                     
              !               -0.47-> N(2D) + O                                                   
              !               -0.03-> N(4S) + O                                                   
              ! -----------                                                                       

              rr = 2.7e-17  ! Cissi Lin [Herron, 1999]                                            

              Reaction = &
                   rr * &
                   Neutrals(iN_2P_) * &
                   Neutrals(iO_3P_)

              IonSources(iNOP_)      = IonSources(iNOP_)      + 0.50 * Reaction
              NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + 0.47 * Reaction
              NeutralSources(iN_4S_) = NeutralSources(iN_4S_) + 0.03 * Reaction
              NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.50 * Reaction

              NeutralLosses(iN_2P_)  = NeutralLosses(iN_2P_)  + Reaction




              ! ----------------------------------------------------------
              ! N(2D)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + O2 -> NO + O(3P) + 3.76 eV
              ! N(2D) + O2 -> NO + O(1D) + 1.80 eV
              ! -----------

              rr = 6.2e-18 *(Tn/300)

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO2_)

              if (UseNeutralConstituent(iO_1D_)) then

                 ! changed (iO_3P_) from 0.9 to 0.98 to reflect Yonker's implementation
                 ! changed (iO_1D_) from 0.1 to 0.02
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + 0.98 * Reaction
                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) + 0.02 * Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.76 * 0.98 + &
                      Reaction * 1.80 * 0.02

                 ChemicalHeatingS(in2d_o2) =  &
                      ChemicalHeatingS(in2d_o2) + &
                      Reaction * 3.76 * 0.98 + &
                      Reaction * 1.80 * 0.02

              else

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 3.76

              endif
              NeutralSources(iNO_)   = NeutralSources(iNO_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iO2_)   = NeutralLosses(iO2_)   + Reaction


              ! Emission caused by Higher Vibrational State from N(2D) + O2 Reaction
              if (UseNOHigherVibrationalStates) then

                  do iVL = 1, 10

                      Emission(i5300_) = Emission(i5300_) + &
                         iVL * Planck_Constant * Speed_Light / 5.3e-6 * Adv1(iVL) / &
                         (Adv1(iVL) + Adv2(iVL) + VSP_f * VSP_k(iVL) * Neutrals(iO_3P_)) &
                         * VSP_g2d(iVL+1) * Reaction

                  enddo

              endif




              ! -----------
              ! N(2D) + NO -> N2 + O + 5.63 eV
              ! -----------

              !rr = 7.0e-17
              rr = 6.7e-17  ! Cissi Lin: Yonker 2013                                              

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iNO_)

              NeutralSources(iN2_)  = NeutralSources(iN2_)  + Reaction
              NeutralSources(iO_3P_)   = NeutralSources(iO_3P_)   + Reaction
              NeutralLosses(iN_2D_) = NeutralLosses(iN_2D_) + Reaction
              NeutralLosses(iNO_)   = NeutralLosses(iNO_)   + Reaction

              ChemicalHeatingSub = &
                   ChemicalHeatingSub + &
                   Reaction * 5.63

              ! ----------------------------------------------------------
              ! O(1D)
              ! ----------------------------------------------------------

              if (UseNeutralConstituent(iO_1D_)) then
                 ! ------------
                 ! O(1D) -> O(3P) + 6300A
                 ! ------------

                 rr = 0.0071

                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)

                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction

                 Emission(iE6300_) = Emission(iE6300_) + Reaction


                 ! ------------
                 ! O(1D) -> O(3P) + 6364A
                 ! ------------

                 rr = 0.0022

                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)

                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction

                 Emission(iE6364_) = Emission(iE6364_) + Reaction

                 ! ------------
                 ! O(1D) + e -> O(3P) + e
                 ! ------------

                 rr = 2.6e-17 * te22m05

                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_)

                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_)  = NeutralLosses(iO_1D_)  + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96

                 ChemicalHeatingS(io1d_e) = &
                      ChemicalHeatingS(io1d_e) + &
                      Reaction * 1.96

                 ! ------------
                 ! O(1D) + N2 -> O(3P) + N2 + 1.96 eV
                 ! ------------

                 rr = 2.3e-17
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iN2_)

                 !              !We create and loose the same amount of N2
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96

                 ChemicalHeatingS(io1d_n2) = &
                      ChemicalHeatingS(io1d_n2) + &
                      Reaction * 1.96


                 ! ------------
                 ! O(1D) + O2 -> O(3P) + O2 + 1.96 eV
                 ! ------------

                 rr = 2.9e-17 * tn06

                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iO2_)

                 ! We create and loose the same amount of O2
                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96

                 ChemicalHeatingS(io1d_o2) = &
                      ChemicalHeatingS(io1d_o2) + &
                      Reaction * 1.96

                 ! ------------
                 ! O(1D) + O(3P) -> O(3P) + O(3P) + 1.96 eV
                 ! ------------

                 rr = 8.0e-18
                 Reaction = &
                      rr * &
                      Neutrals(iO_1D_) * &
                      Neutrals(iO_3P_)

                 NeutralSources(iO_3P_) = NeutralSources(iO_3P_) + Reaction
                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) + Reaction

                 ChemicalHeatingSub = &
                      ChemicalHeatingSub + &
                      Reaction * 1.96

                 ChemicalHeatingS(io1d_o) = &
                      ChemicalHeatingS(io1d_o) + &
                      Reaction * 1.96

              endif

              ! ----------------------------------------------------------
              ! NO
              ! ----------------------------------------------------------
              ! -----------
              ! NO -> NO+ + e
              ! -----------

              !              rr = 6.0e-7

              rr=5.88e-7*(1+0.2*(f107-65)/100)*exp(-2.115e-18* &
                   (Neutrals(iO2_)*1.e-6)**0.8855)

              Reaction = &
                   rr * &
                   Neutrals(iNO_)

              IonSources(iNOP_)   = IonSources(iNOP_)   + Reaction
              NeutralLosses(iNO_) = NeutralLosses(iNO_) + Reaction








!---------------- Production ---------------------

              if (UseNOChemiEquilibrium) then
              ! ----------------------------------------------------------
              ! N2(A)
              ! - added by C Lin 2016
              ! ----------------------------------------------------------

              ! -----------------
              ! What's temporarily missing here is the N2A production rate
              ! -----------------
              ! needs to be added here
              !


              !
              ! Loss of N2(A) - added by C Lin 2016
              !
              ! N2(A) + O(3P) -> NO + N(2D)
              ! rr = 1e-19 * sqrt(T/298) for NOx1d SNOE
              !      [Dilecce and DeBenedictis, 1999]
              ! no T dep. for NOx1D EUVAC [Hill et al., 2000]
              !

              if (UseN2AProfiles) then
                  rr = 1e-19 * sqrt(tn/298.0)
                  Reaction = rr * &
                             Neutrals(iN2_A_) * &
                             Neutrals(iO_3P_)

                  NeutralLosses(iN2_A_)  = NeutralLosses(iN2_A_) + Reaction
                  NeutralLosses(iO_3P_)  = NeutralLosses(iO_3P_) + Reaction
                  NeutralSources(iNO_)   = NeutralSources(iNO_) + Reaction
                  NeutralSources(iN_2D_) = NeutralSources(iN_2D_) + Reaction

                  NOP = NOP + Reaction

              endif



              ! -----------
              ! O2+ + N2 -> NO+ + NO + 0.9333 eV
              ! -----------

              rr = 5.0e-22

              Reaction = &
                   rr * &
                   Ions(iO2P_) * &
                   Neutrals(iN2_)

              NeutralSources(iNO_) = NeutralSources(iNO_) + Reaction

              NOP = NOP + Reaction



              ! N+ + O2 -> O+(4S) + NO + 2.31 eV
              ! -----------

              rr = 3.0e-17

              Reaction = &
                   rr * &
                   Ions(iNP_) * &
                   Neutrals(iO2_)

              NOP = NOP + Reaction



              ! -----------
              ! O+(2P) + N2 -> N+ + NO + 0.70 eV
              ! -----------

              rr = 1.0e-16

              Reaction = &
                   rr * &
                   Ions(iO_2PP_) * &
                   Neutrals(iN2_)

              NOP = NOP + Reaction



              ! -----------
              ! N(4S) + O2 -> NO + O + 1.385 eV
              ! -----------

              !rr = 4.4e-18 * exp(-3220/tn)  ! commented out by C Lin
              rr = 1.5e-17 * exp(-3600/tn)  ! [Sultanou et al., 2006]

              Reaction = &
                   rr * &
                   Neutrals(iN_4S_) * &
                   Neutrals(iO2_)

              NOP = NOP + Reaction






              ! ----------------------------------------------------------
              ! N(2D)
              ! ----------------------------------------------------------

              ! -----------
              ! N(2D) + O2 -> NO + O(3P) + 3.76 eV
              ! N(2D) + O2 -> NO + O(1D) + 1.80 eV
              ! -----------

              rr = 6.2e-18 *(Tn/300)

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_) * &
                   Neutrals(iO2_)

              NOP = NOP + Reaction




!------------Loss--------------------------------------------

              ! N2+ + NO -> NO+ + N2 + 6.33 eV

              rr = 3.3e-16

              Reaction = &
                   rr * &
                   Ions(iN2P_)

              NOL = NOL + Reaction




              ! -----------
              ! O2+ + NO -> NO+ + O2 + 2.813 eV
              ! -----------

              rr = 4.4e-16

              Reaction = &
                   rr * &
                   Ions(iO2P_)

              NOL = NOL + Reaction



              ! -----------
              ! O+(4S) + NO -> NO+ + O + 4.36 eV
              ! -----------

              !!              rr = 8.0e-19
              rr = k3_no_point

              Reaction = &
                   rr * &
                   Ions(iO_4SP_)

              NOL = NOL + Reaction



              ! -----------
              ! NO -> N(4S) + O
              ! -----------

              !              rr = 8.3e-6
              rr=4.5e-6*exp(-1.e-8*(Neutrals(iO2_)*1.e-6)**0.38)

              Reaction = rr

              NOL = NOL + Reaction


              ! -----------
              ! N(4S) + NO -> N2 + O + 3.25 eV
              ! -----------
              !rr = 1.5e-18 * sqrt(tn)  ! commented out by C Lin
              ! rr = 3.4e-17
              rr = 2.1e-17 * exp(100/tn)

              Reaction = &
                   rr * &
                   Neutrals(iN_4S_)

              NOL = NOL + Reaction




              ! -----------
              ! N(2D) + NO -> N2 + O + 5.63 eV
              ! -----------

              rr = 7.0e-17

              Reaction = &
                   rr * &
                   Neutrals(iN_2D_)

              NOL = NOL + Reaction



              ! ----------------------------------------------------------
              ! NO
              ! ----------------------------------------------------------
              ! -----------
              ! NO -> NO+ + e
              ! -----------

              !              rr = 6.0e-7

              rr=5.88e-7*(1+0.2*(f107-65)/100)*exp(-2.115e-18* &
                   (Neutrals(iO2_)*1.e-6)**0.8855)

              Reaction = rr

              NOL = NOL + Reaction

              endif


!-------------------------------------------------------------






              !---- Ions

              if (.not. UseIonChemistry) then
                 IonSources = 0.0
                 IonLosses = 0.0
              else
                 do iIon = 1, nIons-1
                    if (.not.UseIonConstituent(iIon)) then
                       IonSources(iIon) = 0.0
                       IonLosses(iIon) = 0.0
                    endif
                 enddo
              endif

!!!              tli = DtSub * IonLosses
!!!              tsi = DtSub * IonSources + Ions
!!!
!!!              do iIon = 1, nIons-1
!!!                 do while (tsi(iIon)-tli(iIon) < 0.0 .and. DtSub > 1.0e-2)
!!!                    if (tsi(iIon)-tli(iIon) < 0.0 .and. Ions(iIon) < 1.0e7) then
!!!                       IonLosses(iIon) = &
!!!                            (IonSources(iIon) + Ions(iIon)/DtSub)*0.9
!!!                    else
!!!                       DtSub = DtSub/2.0
!!!                    endif
!!!                    tli(iIon) = DtSub * IonLosses(iIon)
!!!                    tsi(iIon) = DtSub * IonSources(iIon) + Ions(iIon)
!!!                 enddo
!!!              enddo
!!!
!!!              !---- Neutrals
!!!
              if (.not. UseNeutralChemistry) then
                 NeutralSources = 0.0
                 NeutralLosses = 0.0
              else
                 do iNeutral = 1, nSpeciesTotal
                    if (.not.UseNeutralConstituent(iNeutral)) then
                       NeutralSources(iNeutral) = 0.0
                       NeutralLosses(iNeutral) = 0.0
                    endif
                 enddo
              endif

              !              if (UseNeutralConstituent(iO_1D_)) then
              !                 NeutralSources(iO_1D_) = NeutralSources(iO_1D_) / y2
              !                 NeutralLosses(iO_1D_) = NeutralLosses(iO_1D_) / y2
              !              endif

!!!              tln = DtSub * NeutralLosses
!!!              tsn = DtSub * NeutralSources + 0.25*Neutrals
!!!              do while (minval(tsn-tln) < 0.0 .and. DtSub > 0.5)
!!!                 DtSub = DtSub/2.0
!!!                 tln = DtSub * NeutralLosses
!!!                 tsn = DtSub * NeutralSources + 0.25*Neutrals
!!!              enddo
!!!
!!!              tli = DtSub * IonLosses
!!!              tsi = DtSub * IonSources + 0.25*Ions
!!!              do while (minval(tsi-tli) < 0.0 .and. DtSub > 0.5)
!!!                 DtSub = DtSub/2.0
!!!                 tli = DtSub * IonLosses
!!!                 tsi = DtSub * IonSources + 0.25*Ions
!!!              enddo
!!!
!!!!!!              Ions(nIons) = 0.0
!!!
!!!              do iIon = 1, nIons-1
!!!
!!!!!!                 if (Ions(iIon) + &
!!!!!!                      (IonSources(iIon) - IonLosses(iIon)) * DtSub < 0.0) then
!!!!!!                    !!!!!! Solve Steady-State !!!!!!!
!!!!!!                    Ions(iIon) = IonSources(iIon)*Ions(iIon)/IonLosses(iIon)
!!!!!!                 else
!!!!!!                    Ions(iIon) = Ions(iIon) + &
!!!!!!                         (IonSources(iIon) - IonLosses(iIon)) * DtSub
!!!!!!                 endif
!!!
!!!!                 Ions(iIon) = max(0.01,Ions(iIon))
!!!
!!!                 if (Ions(iIon) + &
!!!                      (IonSources(iIon) - IonLosses(iIon)) * DtSub < 0.0) then
              ! Take Implicit time step
              Ions(ie_) = 0.0
              do iIon = 1, nIons-1
                 ionso = IonSources(iIon)
                 ionlo = IonLosses(iIon)/Ions(iIon)
                 Ions(iIon) = (Ions(iIon) + ionso * DtSub) / &
                      (1 + DtSub * ionlo)
!!!                 else
!!!                    Ions(iIon) = Ions(iIon) + &
!!!                         (IonSources(iIon) - IonLosses(iIon)) * DtSub
!!!                 endif


                 ! sum for e-
                 Ions(ie_) = Ions(ie_) + Ions(iIon)

                 !                 if (Ions(iIon) < 0.0) then
!!!                 if (isnan(Ions(iIon))) then
!!!                    write(*,*) "Negative Ion Density : ", &
!!!                         iIon, iLon, iLat, iAlt, &
!!!                         Ions(iIon), &
!!!                         IonSources(iIon), IonLosses(iIon)
!!!                 endif
              enddo

              do iNeutral = 1, nSpeciesTotal

!!!                 if (Neutrals(iNeutral) + &
!!!                      (NeutralSources(iNeutral) - &
!!!                       NeutralLosses(iNeutral)) * DtSub < 0.0) then

                 neuso = NeutralSources(iNeutral)
                 neulo = NeutralLosses(iNeutral) / Neutrals(iNeutral)

                 Neutrals(iNeutral)=(Neutrals(iNeutral) + neuso * DtSub) / &
                      (1 + DtSub * neulo)

!!!                 else
!!!
!!!                    Neutrals(iNeutral) = &
!!!                         Neutrals(iNeutral) + &
!!!                         (NeutralSources(iNeutral)-NeutralLosses(iNeutral)) * &
!!!                         DtSub
!!!
!!!                 endif

                 NeutralSourcesTotal(ialt,iNeutral) = &
                      NeutralSourcesTotal(ialt,iNeutral) + &
                      NeutralSources(iNeutral) * DtSub

                 NeutralLossesTotal(ialt,iNeutral) = &
                      NeutralLossesTotal(ialt,iNeutral) + &
                      NeutralLosses(iNeutral) * DtSub

                 !                 if (Neutrals(iNeutral) < 0.0) then
!!!                 if (isnan(Neutrals(iNeutral))) then
!!!                    write(*,*) "Negative Neutral Density : ", &
!!!                         iNeutral, iLon, iLat, iAlt, DtSub, &
!!!                         Neutrals(iNeutral), &
!!!                         NeutralSources(iNeutral), NeutralLosses(iNeutral)
!!!                 endif


              enddo

              ChemicalHeatingRate(iLon,iLat,iAlt) = &
                   ChemicalHeatingRate(iLon,iLat,iAlt) + &
                   ChemicalHeatingSub * DtSub

              ChemicalHeatingSpecies(iLon,iLat,iAlt,:) = &
                   ChemicalHeatingSpecies(iLon,iLat,iAlt,:) + &
                   ChemicalHeatingS * DtSub

              EmissionTotal = EmissionTotal + Emission(:)*DtSub

              DtTotal = DtTotal + DtSub

              if (DtSub < DtMin) DtMin = DtSub

              if (DtSub < 1.0e-9 .and. abs(DtTotal-Dt) > DtSub) then
                 write(*,*) "Chemistry is too fast!!", DtSub

                 ! Check Ions
                 do iIon = 1, nIons
                    write(*,*) "Ion Source/Loss : ", &
                         iIon, IonSources(iIon), IonLosses(iIon)
                 enddo
                 do iNeutral = 1, nSpeciesTotal
                    write(*,*) "Neutral Source/Loss : ", iAlt, &
                         iNeutral, NeutralSources(iNeutral), &
                         NeutralLosses(iNeutral), Neutrals(iNeutral)
                 enddo

                 call stop_gitm("Chemistry is too fast!!")
              endif

              nIters = nIters + 1

           enddo

           IDensityS(iLon,iLat,iAlt,:,iBlock) = Ions
           NDensityS(iLon,iLat,iAlt,:,iBlock) = Neutrals
           if (UseNOChemiEquilibrium) NDensityS(iLon,iLat,iAlt,iNO_,iBlock) = NOP / NOL

           ! modify this for outputing NO emission only - C Lin
           Emissions(iLon, iLat, iAlt, :, iBlock) =  &
                Emissions(iLon, iLat, iAlt, :, iBlock) + EmissionTotal
           !Emissions(iLon, iLat, iAlt, i5300_, iBlock) = NO_Emission

        enddo
     enddo
  enddo

  if (iDebugLevel > 3) then
     do iIon = 1, nIons
        write(*,*) "====> calc_chemistry: Max Ion Density: ", iIon, &
             maxval(IDensityS(1:nLons,1:nLats,(nAlts*4)/5,iIon,iBlock))
     enddo
  endif

  if (iDebugLevel > 2) &
       write(*,*) "===> calc_chemistry: Average Dt for this timestep : ", &
       (Dt*nLats*nLons*nAlts)/nIters

  call end_timing("calc_chemistry_new")

end subroutine calc_chemistry_new
