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

!!! Calculate the expandede grids 
! ------------------------------------------------------------------------
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

!!! Hemispheric integrated values 
! ------------------------------------------------------------------------
subroutine calc_hemi_int_val(nmlt,nmlat,mlat_in,val_in,fac,epr,val_out)
  
  implicit none 

  integer, intent(in) :: nmlt, nmlat
  real, intent(in) :: mlat_in(nmlat+1)
  real, intent(in) :: val_in(nmlt,nmlat), fac, epr
  real, intent(out) :: val_out
  
  real :: colat(nmlat+1), vals(nmlat), cosins(nmlat), val
  integer :: imlat, iimlat
  real, parameter :: pi = 3.14159265359
  
  ! Calculate the average at each MLAT and flip the array
  do imlat=1,nmlat

     iimlat=nmlat-imlat+1
     vals(iimlat) = sum(val_in(:,imlat))/nmlt

  end do

  ! Calculate the colat and flip the array
  colat=90-mlat_in
  colat=colat(nmlat+1:1:-1)*epr

  ! Calculate the cosin difference
  colat=colat/180.*pi
  cosins=cos(colat(1:nmlat))-cos(colat(2:nmlat+1))

  ! Caculate the sum of each latitude
  val_out=0.
      
  do imlat=1,nmlat

     val=vals(imlat)*cosins(imlat)*2*pi
     if (val<0) val=1.0e-20
         
     val_out=val_out+val
  end do
      
  ! Multiply the factor
  val_out=val_out*(6500**2)*fac

end subroutine calc_hemi_int_val

!!! Hemispheric integrated diff_ef
! ------------------------------------------------------------------------
subroutine calc_hemi_int_val_channel(nmlt,nmlat,nchannel,diff_ef,mlat_in,&
     fac,epr,val_out)
  
  implicit none
  
  integer, intent(in) :: nmlt, nmlat, nchannel
  real, intent(in) :: diff_ef(nmlt,nmlat,nchannel)
  real, intent(in) :: mlat_in(nmlat+1), fac, epr
  real, intent(out) :: val_out(nchannel)
  
  integer :: ichannel
  
  do ichannel=1,nchannel
     
     call calc_hemi_int_val(nmlt,nmlat,mlat_in,&
          diff_ef(:,:,ichannel),fac,epr,val_out(ichannel))
     
  end do

end subroutine calc_hemi_int_val_channel
      
!!! Calculate the total energy and number fluxes 
! ------------------------------------------------------------------------
subroutine calc_total_ef_nf(nmlt,nmlat,nchannel,channels,diff_ef,&
     total_ef,total_nf)
  
  implicit none 
  
  integer, intent(in) :: nmlt, nmlat, nchannel
  real, intent(in) :: channels(nchannel)
  real, intent(in) :: diff_ef(nmlt,nmlat,nchannel)
  real, intent(out) :: total_ef(nmlt,nmlat), total_nf(nmlt,nmlat)
  
  real :: vals(nchannel)
  integer :: imlt, imlat
  real, parameter :: pi = 3.14159265359

  do imlt=1,nmlt
     do imlat=1,nmlat

        ! qingyu, 10/23/2019
        ! Use Channels >500 eV for the calculation
        ! Channels 1:11
        ! According to Robinson et al., (1987) and our comparisons
        
        vals=diff_ef(imlt,imlat,:)

        total_nf(imlt,imlat) = sum(vals(1:11))*(1.0e8)*pi*0.364 ! cm-2 s-1
        total_ef(imlt,imlat) = dot_product(vals(1:11),channels(1:11)) * &
             (1.0e-4)*pi*0.364*1.602 ! erg cm-2 s-1

     end do
  end do

end subroutine calc_total_ef_nf

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

