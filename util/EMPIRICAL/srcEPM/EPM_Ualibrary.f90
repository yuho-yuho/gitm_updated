! Higher-class functions for the EPM
! Updated from ModSEP_Ualibrary.f90
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!! Reconstruct the SPH coefficients 
! -----------------------------------------------------------------------------
subroutine reconst_spharm_coeffs(coeff_in,ca_in,LMAX,NMAX,coeff_out)

  implicit none 

  integer, intent(in) :: LMAX, NMAX
  real, intent(in) :: coeff_in((LMAX+1)**2,2*NMAX+1) 
  real, intent(in) :: ca_in
  real, intent(out) :: coeff_out((LMAX+1)**2)

  integer :: iparm                                      
  real :: PARM(2*NMAX+1), fs(2*NMAX+1)

   do iparm=1,(LMAX+1)**2                                                    
      PARM(:) = coeff_in(iparm,:)                                            
      call construct_fourier_series(ca_in,NMAX,fs)                           
      coeff_out(iparm)=dot_product(PARM,fs)                                  
   end do

end subroutine reconst_spharm_coeffs

!!! Reconstruct the electric potential at a certain point 
! -----------------------------------------------------------------------------
subroutine reconst_epot_sp(coeff_in,LMAX,theta,phi,epot)

  use modspharm, only: reconst_spharm
  use epm_interface, only: pi

  implicit none                                                             
  
  integer, intent(in) :: LMAX
  real, intent(in) :: coeff_in((LMAX+1)**2)                                 
  real, intent(in) :: theta, phi                                            
  real, intent(out) :: epot                                                 
  
  real :: spharms((LMAX+1)**2), theta1

  ! This is to remove the sharp boundary at 45 MLAT
  if (theta>pi) then
     theta1=pi
  else
     theta1=theta
  end if

  call reconst_spharm(theta1,phi,LMAX,spharms)
  epot = dot_product(coeff_in,spharms)

end subroutine reconst_epot_sp

!!! Obtain the slope and yint for a given ca_in
subroutine reconst_epm_slope_yint(coeff1_in,coeff2_in,NMAX,ca_in,slope,yint)
  
  use EPM_Interface

  implicit none 

  integer, intent(in) :: NMAX
  real, intent(in) :: coeff1_in(2*NMAX+1), coeff2_in(2*NMAX+1), ca_in
  real, intent(out) :: slope, yint

  real :: fs(2*NMAX+1)

  call construct_fourier_series(ca_in,NMAX,fs)

  if (epm_debug) write(*,*) "Slope Coefficients:", coeff1_in

  slope=dot_product(coeff1_in,fs)
  yint=dot_product(coeff2_in,fs)

  slope=slope*(1e-4)

end subroutine reconst_epm_slope_yint

!!! Correct the CPCP
subroutine correct_cpcp(cf_in,cf_inf,slope,yint,cpcp_in,scale)

  use epm_interface, only: epm_debug

  implicit none 

  real, intent(in) :: cf_in, cf_inf, slope, yint, cpcp_in
  real, intent(out) :: scale

  real :: beta_in, cpcp

  call calc_beta(cf_in,cf_inf,beta_in)

  cpcp = slope*beta_in+yint

  if (cpcp_in<=0.) then 
     scale=1.
  else
     scale = cpcp/cpcp_in
  end if

  if (epm_debug) then 
       write(*,*) "-----------------------------------------------"
       write(*,*) slope, yint, cpcp, cpcp_in
       write(*,*) "-----------------------------------------------"
    end if

end subroutine correct_cpcp

!!! Calculate the expansion 
! -----------------------------------------------------------------------------
subroutine calc_epm_expansion(cf_in,cf_ref,ca_in,cf_inf,expansion_rate)

  implicit none 

  real, intent(in) :: cf_in, cf_ref, ca_in, cf_inf
  real, intent(out) :: expansion_rate      
                                                                              
  real :: beta_in, beta_ref                                                   
  real :: half_theta, sin_ht, ca_factor                                       
  real, parameter :: pi = 3.14159265359

  call calc_beta(cf_in,cf_inf,beta_in)                                         
  call calc_beta(cf_ref,cf_inf,beta_ref)

  half_theta=(ca_in)/2.                                                        
  sin_ht=sin(half_theta)                                                       
  ca_factor=sin_ht**2

  expansion_rate=1+ca_factor*&
       ((9.31*(1e-4)*beta_in+21.4)/(9.31*(1e-4)*beta_ref+21.4)-1)

end subroutine calc_epm_expansion

!!! Calculate the displacement 
! -----------------------------------------------------------------------------
subroutine calc_epm_displacement(cf_in,cf_inf,ca_in,cf_ref,disp)

  real, intent(in) :: cf_in, cf_ref, ca_in, cf_inf
  real, intent(out) :: disp      
                                                                              
  real :: beta_in, beta_ref                                                   
  real :: theta, sin_t, ca_factor                                       
  real, parameter :: pi = 3.14159265359, disp0=1.0e-4

  call calc_beta(cf_in,cf_inf,beta_in)                                         
  call calc_beta(cf_ref,cf_inf,beta_ref)

  theta=ca_in
  sin_t=sin(theta)

  disp=disp0*sin_t*(beta_in-beta_ref)

  if (disp>1.5) disp=1.5
  if (disp<-1.5) disp=-1.5

end subroutine calc_epm_displacement

!!! calculate the theta and phi for a given mlat and mlt 
! -----------------------------------------------------------------------------
subroutine calc_theta_phi(mlat_in,mlt_in,er,disp,theta,phi)
  
  use epm_interface, only: pi

  implicit none                                                             
  
  real, intent(in) :: mlat_in, mlt_in, er, disp                            
  real, intent(out) :: theta, phi                                          
  
  real :: mlt1, mlat1                                                      
  real :: t1, r1, t2, r2
  real :: x1, y1, x2, y2

  mlt1=mlt_in
  if (mlt_in<0.) mlt1=mlt_in+24.                                            
  if (mlt_in>24.) mlt1=mlt_in-24.

  mlat1=abs(mlat_in)

  t1=(mlt1-6)/12*pi                                                         
  r1=(90-mlat1)                                                             
  
  x1=r1*cos(t1)                                                             
  y1=r1*sin(t1)

  x2=x1/er-disp                                                             
  y2=y1/er

  r2=sqrt(x2**2+y2**2)                                                      
  t2=atan2(y2,x2)

  t2=t2+pi/2.                                                               
  
  if (t2<0.) t2=t2+2*pi                                                     
  
  phi=t2                                                                    
  theta=r2/45.*pi

end subroutine calc_theta_phi

!!! Get the coefficients, extrapolation, expansion and displacement 
! -----------------------------------------------------------------------------
subroutine calc_epm_coeffs(cf_in,cf_inf,cf_inf1,ncf,ref_cfs,cf0,ca_in,&
     LMAX,NMAX,coeff,er,disp)

  use epm_interface

  implicit none 

  integer, intent(in) :: ncf, LMAX, NMAX
  real, intent(in) :: cf_in, cf_inf, cf_inf1, ca_in, ref_cfs(ncf),cf0
  real, intent(out) :: coeff((LMAX+1)**2), er, disp

  integer :: imod1, imod2                                                   
  real :: wgt1, wgt2                                                        
  real :: coeff1((LMAX+1)**2), coeff2((LMAX+1)**2)

  !! Determine coeff
  call determine_cf_models(cf_in,ncf,ref_cfs,imod1,imod2)
  call determine_epm_weights(cf_in,cf_inf,cf0,ncf,ref_cfs,&
       imod1,imod2,wgt1,wgt2)

  if (epm_debug) write(*,*) imod1, imod2, wgt1, wgt2

  ! Coeff1
  if (imod2==1) then
     coeff1(:)=epm_cat0_coeffs(:)
  else
     call reconst_spharm_coeffs(epm_all_ca_coeffs(imod1,:,:),ca_in,&
          LMAX,NMAX,coeff1)
  end if

  ! Coeff2
  if (imod2==ncf) then
     coeff2(:) = coeff1(:)
  else
     call reconst_spharm_coeffs(epm_all_ca_coeffs(imod2,:,:),ca_in,&
          LMAX,NMAX,coeff2)
  end if

  ! Combine coeff1 and coeff2 --> coeff
  coeff(:) = coeff1(:) * wgt1 + coeff2(:) * wgt2 
  
  !! Determine er and disp
  er=1.
  disp=0.

  if (imod1==ncf) then
     call calc_epm_expansion(cf_in,ref_cfs(ncf),ca_in,cf_inf1,er)
     call calc_epm_displacement(cf_in,cf_inf,ca_in,ref_cfs(ncf),disp)
  end if

  if (epm_debug) write(*,*) "------ er, disp ------", er, disp

end subroutine calc_epm_coeffs

!!! Calculate the electric potential 
! -----------------------------------------------------------------------------
subroutine calc_epm_epot(cf_in,cf_inf,cf_inf1,ca_in,LMAX,NMAX,nmlt,mlts,&
     nmlat,mlats,coeff,scale,er,disp)

  use epm_interface 

  implicit none 

  real, intent(in) :: cf_in,cf_inf,cf_inf1,ca_in
  integer, intent(in) :: nmlt, nmlat, LMAX, NMAX
  real, intent(in) :: mlts(nmlt), mlats(nmlat)
  real, intent(out) :: scale, er, disp
  real, intent(out) :: coeff((LMAX+1)**2)

  integer :: imlt, imlat
  real :: epot(nmlt,nmlat)
  real :: mlat1, mlt1, theta, phi, pot, slope, yint, cpcp1, cpcp2
  real :: maxi, mini

  ! Calculate coeffcients 
  call calc_epm_coeffs(cf_in,cf_inf,cf_inf1,epm_ncf,epm_ref_cfs,epm_cf0,ca_in,&
       LMAX,NMAX,coeff,er,disp) 

  ! Calculate gridded electric potential 
  do imlt=1,nmlt
     
     mlt1=mlts(imlt)

     do imlat=1,nmlat
        
        mlat1=mlats(imlat)
        call calc_theta_phi(mlat1,mlt1,er,disp,theta,phi) 
        call reconst_epot_sp(coeff,LMAX,theta,phi,pot)
        epot(imlt,imlat)=pot

        !if (epm_debug) write(*,*) mlt1, mlat1, theta, phi, pot

     end do

  end do

  ! Calculate the cpcp of the epot
  call find_max_min(epot,nmlt,nmlat,maxi,mini)
  cpcp1=maxi-mini
  if (epm_debug) write(*,*) maxi,mini,"CPCP1 (kV):", cpcp1

  ! Calculate slope and yint

  if (epm_debug) then 

     write(*,*) "Slope coeffs", epm_slope_coeffs

  end if

  call reconst_epm_slope_yint(epm_slope_coeffs,epm_yint_coeffs,4,ca_in,&
       slope,yint) 

  ! Calculate the correction 
  call correct_cpcp(cf_in,cf_inf,slope,yint,cpcp1,scale)
  if (epm_debug) write(*,*) "-------- Scale --------", scale

end subroutine calc_epm_epot

!!! Caculate the electric potential at a certain location 
! -----------------------------------------------------------------------------
subroutine calc_epm_sp_epot(mlt,mlat,coeff,scale,er,disp,pot)

  use epm_interface, only: epm_LMAX
  
  implicit none 

  real, intent(in) :: mlt, mlat, scale, er, disp
  real, intent(in) :: coeff((epm_LMAX+1)**2)
  real, intent(out) :: pot

  real :: theta, phi

  call calc_theta_phi(mlat,mlt,er,disp,theta,phi)                       
  call reconst_epot_sp(coeff,epm_LMAX,theta,phi,pot)

  pot = pot * scale

end subroutine calc_epm_sp_epot
