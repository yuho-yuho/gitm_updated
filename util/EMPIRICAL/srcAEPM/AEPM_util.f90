! Provide basic functions for the AEPM 
! Based on the ModSAM_Ualibrary.f90
! Created: Qingyu Zhu, 05/17/2020
!
! -----------------------------------------------------------------------------

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
subroutine determine_weights(cf_in,cf0,ncf,ref_cfs,imod1,imod2,wgt1,wgt2)

  implicit none 

  real, intent(in) :: cf_in, cf0
  integer, intent(in) :: ncf, imod1, imod2 
  real, intent(in) :: ref_cfs(ncf)
  real, intent(out) :: wgt1, wgt2

  real :: cf1, cf2
  
  ! Small cf_in
  if (imod2==1) then

     if (cf_in<cf0) then
        wgt1=1.
        wgt2=0.
     else
        cf1=cf0
        cf2=ref_cfs(imod2)
        
        wgt1=(cf2-cf_in)/(cf2-cf1)
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
     
     wgt1=(cf2-cf_in)/(cf2-cf1)
     wgt2=1-wgt1

  end if

end subroutine determine_weights

