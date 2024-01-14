! This module provides very basic functions for the EPM module
! Created: Qingyu Zhu, 05/27/2020
!
! ----------------------------------------------------------------------------

!!! Removed some functions that can be found in 
!!! AEPM_util_base.f90 or AEPM_util.f90
!!! Determine the weight of the models used for the final pattern
! -----------------------------------------------------------------------------
subroutine determine_epm_weights(cf_in,cf_inf,cf0,ncf,ref_cfs,&
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

end subroutine determine_epm_weights

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
