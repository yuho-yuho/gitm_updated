! This module provides very basic functions for the AEPM module
! Created: Qingyu Zhu, 05/18/2020
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

!!! Calculate the expanded MLTs
subroutine calc_expanded_grid(nmlat,mlat_in,epr,mlat_out)                      
                                                                               
  implicit none                                                                
                                                                               
  integer, intent(in) :: nmlat                                                 
  real, intent(in) :: mlat_in(nmlat), epr                                      
  real, intent(out) :: mlat_out(nmlat)                                         
                                                                               
  integer :: ii                                                                
  real :: colat                                                                
                                                                               
  do ii=1, nmlat                                                               
                                                                               
     colat=(90-mlat_in(ii))*epr                                                
     mlat_out(ii)=90.-colat                                                    
                                                                               
  end do                                                                       
                                                                               
end subroutine calc_expanded_grid

!!! Find the the location of a GITM MLT grid 
! ------------------------------------------------------------------------
subroutine find_mlt_loc(mlt_in,mlts,nmlt,loc1,loc2,wgt1,wgt2)
  
  implicit none
  
  integer, intent(in) :: nmlt
  real, intent(in) :: mlt_in, mlts(nmlt)
  integer, intent(out) :: loc1, loc2
  real, intent(out) :: wgt1, wgt2
  
  real :: mlt1, mlt2
  integer :: ii

  loc1=-99
  loc2=-99
  
  mlt1=0.
  mlt2=1.

  if ((mlt_in<mlts(1)) .or. (mlt_in>=mlts(nmlt)) ) then
     
     loc1 = nmlt
     loc2 = 1

     mlt1 = mlts(loc1)
     mlt2 = mlts(loc2)

     if ((mlt_in-mlt2)<0) then
        mlt1=mlt1-24
     else
        mlt2=mlt2+24
     end if

  else

     do ii=1,nmlt-1

        if ((mlt_in>=mlts(ii)) .and. (mlt_in<mlts(ii+1))) then
           loc1=ii
           loc2=ii+1
           mlt1=mlts(loc1)
           mlt2=mlts(loc2)
           exit
        end if

     end do

  end if

  wgt1=(mlt2-mlt_in)/(mlt2-mlt1)
  wgt2=1-wgt1
  
end subroutine find_mlt_loc

!!! Find the location of a GITM MLAT grid
! -----------------------------------------------------------------------------
subroutine find_mlat_loc(mlat_in,mlats,nmlat,loc1,loc2,wgt1,wgt2)

  implicit none 

  integer, intent(in) :: nmlat
  real, intent(in) :: mlat_in, mlats(nmlat)
  integer, intent(out) :: loc1, loc2
  real, intent(out) :: wgt1, wgt2
  
  real :: mlat1, mlat2
  integer :: ii
      
  loc1=-99
  loc2=-99

  mlat1=50.
  mlat2=51.

  wgt1=1.
  wgt2=0.

  if (mlat_in<mlats(1)) then
     
     loc1=1
     loc2=1
     wgt1=0.
     wgt2=1.

  else if (mlat_in>=mlats(nmlat)) then

     loc1=nmlat
     loc2=nmlat
     wgt1=1.
     wgt2=0.
         
  else

     do ii=1,nmlat-1
        if ((mlat_in>=mlats(ii)) .and. (mlat_in<mlats(ii+1))) then
           loc1=ii
           loc2=ii+1
           mlat1=mlats(loc1)
           mlat2=mlats(loc2)
           exit
        end if
     end do

     wgt1=(mlat2-mlat_in)/(mlat2-mlat1)
     wgt2=1-wgt1
     
  end if
         
end subroutine find_mlat_loc

!!! Get the value on a GITM grid
! ------------------------------------------------------------------------
subroutine get_grid_value(mlt_in,mlat_in,mlts,mlats,nmlt,nmlat,&
     val_in,val_out)
  
  implicit none 
  
  real, intent(in) :: mlt_in, mlat_in
  integer, intent(in) :: nmlt, nmlat
  real, intent(in) :: mlts(nmlt), mlats(nmlat), val_in(nmlt,nmlat)
  real, intent(out) :: val_out
  
  integer :: imlt1, imlt2, imlat1, imlat2
  real :: wgt11, wgt12, wgt21, wgt22
  real :: val11, val12, val21, val22
  
  call find_mlt_loc(mlt_in,mlts,nmlt,imlt1,imlt2,wgt11,wgt12)
  call find_mlat_loc(mlat_in,mlats,nmlat,imlat1,imlat2,wgt21,wgt22)
  
  val11 = val_in(imlt1,imlat1)*wgt11*wgt21
  val12 = val_in(imlt1,imlat2)*wgt11*wgt22
  val21 = val_in(imlt2,imlat1)*wgt12*wgt21
  val22 = val_in(imlt2,imlat2)*wgt12*wgt22
  
  val_out = val11 + val12 + val21 + val22
      
end subroutine get_grid_value

