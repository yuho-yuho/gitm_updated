! Following subroutines provide the new calculations of 
! the collision frequencies
!
! History:
!    Created: 07/05/17, Qingyu Zhu, qingyu.zhu@mavs.uta.edu
!
! Note:
!    Formulas come from Richmond (2014), Ionospheric Electrodynamics, Space
! Weather fundamentals, Ch 14
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine calc_new_collisions(iBlock)

! This subroutine provides the calculation for different ion-neutral
! collision frequencies

  use ModGITM
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock

  real, parameter :: NOP_N2 = 4.35, NOP_O2 = 4.35, NOP_O = 1.9,&
       O2P_N2 = 4.3, O2P_O2 = 5.2, O2P_O = 1.8,&
       OP_N2 = 5.4, OP_O2 = 7.0, OP_O = 6.7


  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       Tn, Ti, Tr, Tr1, Nden_N2, Nden_O2, Nden_O, &
       IRho_NOP, IRho_O2P, IRho_OP


  ! calculate Ti, Tn, Tr, deal with the corners
  
  Tn = Temperature(:,:,:,iBlock) * TempUnit(:,:,:)
  Ti = ITemperature(:,:,:,iBlock)
  Tr = 0.5*(Tn+Ti)

  Tr1 = Tr
  where(Tr1<235.) Tr1 = 235.

  Tr1 = Tr1/1000.

  Nden_N2 = NDensityS(:,:,:,iN2_,iBlock)
  Nden_O2 = NDensityS(:,:,:,iO2_,iBlock)
  Nden_O = NDensityS(:,:,:,iO_3P_,iBlock)

  IonCollisions1(:,:,:,INOP_) = 1e-16*((NOP_N2*Nden_N2+&
       NOP_O2*Nden_O2)/Tr1**(0.11)+NOP_O*Nden_O/Tr1**(0.19))

  IonCollisions1(:,:,:,IO2P_) = 1e-16*(O2P_N2*Nden_N2+&
       O2P_O2*Nden_O2+O2P_O*Nden_O/Tr1**(0.19))

  IonCollisions1(:,:,:,IO_4SP_) = 1e-16*(OP_N2*Nden_N2/Tr1**(-0.20)+&
       OP_O2*Nden_O2*Tr1**(0.05)+&
       OP_O*Nden_O*Tr1**(0.5)*(0.96-0.135*log10(Tr1))**2) 

  ! To compare with previous one, we need to calculate 
  ! the effective collision frequency 

  ! sum(mi*ni*vin)/sum(mi*ni)

  IRho_NOP = MassI(iNOP_) * IDensityS(:,:,:,iNOP_,iBlock)
  IRho_O2P = MassI(iO2P_) * IDensityS(:,:,:,iO2P_,iBlock)
  IRho_OP = MassI(iO_4SP_) * IDensityS(:,:,:,iO_4SP_,iBlock)

  Collisions1(:,:,:,iVIN_) = 0.
  Collisions1(:,:,:,iVIN_) = (IRho_NOP*IonCollisions1(:,:,:,INOP_)+&
       IRho_O2P*IonCollisions1(:,:,:,IO2P_)+&
       IRho_OP*IonCollisions1(:,:,:,IO_4SP_))/&
       (IRho_NOP+IRho_O2P+IRho_OP)

end subroutine calc_new_collisions
