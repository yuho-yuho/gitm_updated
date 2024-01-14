! This module provides very basic functions for the EPM module
! Created: Qingyu Zhu, 05/27/2020
!
! ----------------------------------------------------------------------------

!!! Get the coupling function and IMF inputs 
! -------------------------------------------------------------------------
subroutine calc_cf_ca(by,bz,vx,n,cf,ca)                                        
                                                                               
  implicit none                                                                
                                                                               
  real, intent(in) :: by, bz, vx, n                                            
  real, intent(out) :: cf, ca                                                  
                                                                               
  real :: bt  
  real, parameter :: pi = 3.14159265359
                                                                               
  bt=sqrt(by**2+bz**2)                                                         
  cf=((vx)**(4./3.))*((bt)**(2./3.))*((n)**(1./6.))                            
  ca=atan2(by,bz)                                                              
                                                                               
  if (ca<0.) ca=ca+2*pi                                                        
                                                                              
end subroutine calc_cf_ca 

!!! Reconstruct the sines and cosines
! ------------------------------------------------------------------------
subroutine construct_fourier_series(x,order,coeff)
  
  implicit none
  
  real, intent(in) :: x
  integer, intent(in) :: order 
  real, intent(out) :: coeff(2*order+1)
  
  integer :: ii
  real :: t
      
  coeff(1)=1.

  do ii=1,order
     t=ii*x
     coeff(2*ii)=cos(t)
     coeff(2*ii+1)=sin(t)
  end do
  
end subroutine construct_fourier_series

!!! Calculate beta 
! ------------------------------------------------------------------------
subroutine calc_beta(cf_in,cf_inf,beta)

  implicit none 

  real, intent(in) :: cf_in, cf_inf
  real, intent(out) :: beta

  beta=cf_in/sqrt(1+(cf_in/cf_inf)**2)

end subroutine calc_beta

!!! Determine the module will be used for the IMF and solar wind inputs 
! -----------------------------------------------------------------------------
subroutine determine_cf_models(cf_in,ncf,ref_cfs,imod1,imod2)

  implicit none 

  real, intent(in) :: cf_in
  integer, intent(in) :: ncf
  real, intent(in) :: ref_cfs(ncf)
  integer, intent(out) :: imod1, imod2
  
  integer :: ii 

  if (cf_in<ref_cfs(1)) then 
     imod1=1
     imod2=1
  else if (cf_in>=ref_cfs(ncf)) then
     imod1=ncf
     imod2=ncf
  else

     do ii=1,ncf-1
        if ((cf_in>=ref_cfs(ii)) .and. (cf_in<ref_cfs(ii+1))) then
           imod1=ii
           imod2=ii+1
           exit
        end if
     end do

  end if

end subroutine determine_cf_models

!!! Determine the weight of the models used for the final pattern
! -----------------------------------------------------------------------------
subroutine determine_weights(cf_in,cf_inf,cf0,ncf,ref_cfs,&
     imod1,imod2,wgt1,wgt2)

  implicit none 

  real, intent(in) :: cf_in, cf_inf, cf0
  integer, intent(in) :: ncf, imod1, imod2 
  real, intent(in) :: ref_cfs(ncf)
  real, intent(out) :: wgt1, wgt2

  real :: cf1, cf2
  real :: beta_in, beta1, beta2
  
  call calc_beta(cf_in,cf_inf,beta_in)

  ! Small cf_in
  if (imod2==1) then

     if (cf_in<cf0) then
        wgt1=1.
        wgt2=0.
     else
        cf1=cf0
        cf2=ref_cfs(imod2)
        
        call calc_beta(cf1,cf_inf,beta1)
        call calc_beta(cf2,cf_inf,beta2)

        wgt1=(beta2-beta_in)/(beta2-beta1)
        wgt2=1-wgt1
     end if

  ! Large cf_in
  else if (imod1==ncf) then
     wgt1=0.
     wgt2=1.

  ! Intermedium
  else
     cf1=ref_cfs(imod1)
     cf2=ref_cfs(imod2)
     
     call calc_beta(cf1,cf_inf,beta1)                                        
     call calc_beta(cf2,cf_inf,beta2)

     wgt1=(beta2-beta_in)/(beta2-beta1)
     wgt2=1-wgt1

  end if

end subroutine determine_weights

!!! Find the maximum and minimum of an 2d array
! -----------------------------------------------------------------------------
subroutine find_max_min(arr_in,nmlt,nmlat,maxi,mini)

  integer, intent(in) :: nmlt, nmlat
  real, intent(in) :: arr_in(nmlt,nmlat)
  real, intent(out) :: maxi, mini

  real :: arr(nmlt)

  arr=maxval(arr_in,dim=2)
  maxi=maxval(arr)

  arr=minval(arr_in,dim=2)
  mini=minval(arr)

end subroutine find_max_min
