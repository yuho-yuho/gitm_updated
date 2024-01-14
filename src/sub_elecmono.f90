! This is the fortran version of Fang et al. (2010) parameterization method
! to the resulting ionization altitude profile from isotropically 
! precipitating monoenergetic electrons in the Earth's atmosphere.
!
! The original IDL version is Copyrighted by: 
! Xiaohua Fang 
! Laboratory for Atmospheric and Space Physics University of Colorado 
! 392 UCB, Boulder, CO 80309-0392  
! xiaohua.fang@lasp.colorado.edu
!
! Reference:
! Fang, X., C. Randall, D. Lummerzheim, W. Wang, G. Lu, S. Solomon, 
! and R. Frahm (2010) 
! Parameterization of monoenergetic electron impact ionization,
! Geophys. Res. Lett., 37, L22106, doi:10.1029/2010GL045406.
!
! Rewriten in Fortran by:
! Qingyu Zhu, 08/09/2020
!
! -----------------------------------------------------------------------------
subroutine sub_elecmono(Emono, Qmono, rho, H, qtot)

  ! Input(s):
  ! Emono: monoenergetic electron energy [keV]
  !        (valid energy range: 0.1 < Emono < 1.E3  keV)
  ! Qmono: total incident energy flux from monoenergetic electrons 
  !        [erg cm-2 s-1]
  ! rho: altitude profile of Earth's atmospheric mass density [g cm-3]
  ! H: altitude profile of Earth's atmospheric scale height [km]
  
  ! Output(s):
  ! qtot: altitude profile of total resulting ionization rate [cm-3 s-1]

  implicit none 

  ! Dummy arguments
  real, intent(in) :: Emono, Qmono, rho, H 
  real, intent(out) :: qtot

  ! Local variables
  real :: ScaleHeight, IncidentEnergyFlux, y

  real, parameter :: P(4,8) = &
       RESHAPE((/1.24616E+00,  1.45903E+00, -2.42269E-01,  5.95459E-02, &
       2.23976E+00, -4.22918E-07,  1.36458E-02,  2.53332E-03, &
       1.41754E+00,  1.44597E-01,  1.70433E-02,  6.39717E-04, &
       2.48775E-01, -1.50890E-01,  6.30894E-09,  1.23707E-03, &
       -4.65119E-01, -1.05081E-01, -8.95701E-02,  1.22450E-02, &
       3.86019E-01,  1.75430E-03, -7.42960E-04,  4.60881E-04, &
       -6.45454E-01,  8.49555E-04, -4.28581E-02, -2.99302E-03, &
       9.48930E-01,  1.97385E-01, -2.50660E-03, -2.06938E-03/),&
       (/4,8/))

  real :: LnEmono, C(8), LnY, fy
  integer :: i, j

  !! Main function
  ScaleHeight = H*1.e5 ! [km-->cm]

  IncidentEnergyFlux = Qmono*6.2415e8 ! [erg cm-2 s-1 --> keV cm-2 s-1]

  ! Eq. (1)
  y = 2./Emono*(rho*ScaleHeight/6.E-6)**0.7 !Normalized atmospheric column mass
  
  lnEmono = alog(Emono)
  C=0.

  ! Eq. (5)
  do i=1,8
     do j=1,4
        
        C(i)=C(i) + P(j,i) * (LnEmono**(j-1))

     end do
  end do
        
  C=EXP(C)

  ! Eq. (4)
  LnY = alog(y)
  fy = C(1)*EXP( C(2)*lnY - C(3)*EXP(C(4)*lnY) ) + &
       C(5)*EXP( C(6)*lnY - C(7)*EXP(C(8)*lnY) )


  if ((Emono>=0.1) .and. (Emono<=1.e3)) then
     qtot = IncidentEnergyFlux/3.5E-2/ScaleHeight*fy
  else
     qtot =0.
  end if  

end subroutine sub_elecmono
