! Higher-class functions of the EFVM

!! RECONSTRUCT THE MLT COEFFICIENTS
! -----------------------------------------------------------------------------
subroutine reconstruct_efvm_mlt_coeffs(coeff_in,nmlat,LMAX,NMAX,ca_in,&
     coeff_out)

  implicit none

  integer, intent(in) :: nmlat,LMAX,NMAX                              
  real, intent(in) :: coeff_in(nmlat,2*LMAX+1,2*NMAX+1)               
  real, intent(in) :: ca_in                                                    
  real, intent(out) :: coeff_out(nmlat,2*LMAX+1)

  integer :: imlat, iparm
  real :: PARM(2*NMAX+1), fs(2*NMAX+1)

  do imlat=1,nmlat                                                           
     do iparm=1,2*LMAX+1                                                     
        
        PARM(:)=coeff_in(imlat,iparm,:)                             
        call construct_fourier_series(ca_in,NMAX,fs)                         
        coeff_out(imlat,iparm)=dot_product(PARM,fs)                 
        
     end do
  end do

end subroutine reconstruct_efvm_mlt_coeffs

!! RECONSTRUCT THE DISTRIBUTIONS OF THE DED1 AND DED2
! -----------------------------------------------------------------------------
subroutine reconst_efvm_dEd(coeff_in,nmlt,mlts,nmlat,LMAX,dEd_out)

  implicit none 

  integer, intent(in) :: nmlat,LMAX,nmlt                              
  real, intent(in) :: mlts(nmlt)                                               
  real, intent(in) :: coeff_in(nmlat,2*LMAX+1)                        
  real, intent(out) :: dEd_out(nmlt,nmlat)                        
                                                                               
  integer :: imlat, imlt                                             
  real :: PARM(2*LMAX+1), fs(2*LMAX+1), phi, val                               
                                                                               
  real, parameter :: pi = 3.14159265359

  do imlat=1,nmlat

     PARM=coeff_in(imlat,:)

     do imlt=1,nmlt

        phi=mlts(imlt)/12*pi
        call construct_fourier_series(phi,LMAX,fs)
        val=dot_product(PARM,fs)
        
        if (val<0) val=0.
        dEd_out(imlt,imlat)=val

     end do ! MLT

  end do ! MLAT

end subroutine reconst_efvm_dEd

!! CALCULATE THE COEFFICETS
! -----------------------------------------------------------------------------
subroutine calc_efvm_coeffs(cf_in,cf_inf,ncf,ref_cfs,cf0,ca_in,nmlat,&
     LMAX,NMAX,coeff1,coeff2)

  use efvm_interface

  implicit none

  integer, intent(in) :: ncf, LMAX, NMAX, nmlat
  real, intent(in) :: cf_in, cf_inf, ca_in, ref_cfs(ncf),cf0
  real, intent(out) :: coeff1(nmlat,(LMAX+1)**2), coeff2(nmlat,(LMAX+1)**2)

  integer :: imod1, imod2
  real :: wgt1, wgt2

  real :: coeff3(nmlat,(LMAX+1)**2), coeff4(nmlat,(LMAX+1)**2)
  
  ! Determine coeff                                                            
  call determine_cf_models(cf_in,ncf,ref_cfs,imod1,imod2)                      
  call determine_efvm_weights(cf_in,cf_inf,cf0,ncf,ref_cfs,&                
       imod1,imod2,wgt1,wgt2)                                          
                                                                            
  if (efvm_debug) write(*,*) imod1, imod2, wgt1, wgt2 

  !! Ed1
  ! Coeff3
  if (imod2==1) then
     coeff3(:,:)=efvm_cat0_coeffs1(:,:)
  else
     call reconstruct_efvm_mlt_coeffs(efvm_all_ca_coeffs1(imod1,:,:,:),nmlat,&
          LMAX,NMAX,ca_in,coeff3)
  end if

  ! Coeff4
  if (imod2==ncf) then
     coeff4(:,:)=coeff3(:,:)
  else
     call reconstruct_efvm_mlt_coeffs(efvm_all_ca_coeffs1(imod2,:,:,:),nmlat,& 
          LMAX,NMAX,ca_in,coeff4)                                              
  end if

  coeff1(:,:) = coeff3(:,:) * wgt1 + coeff4(:,:) * wgt2 

  !! Ed2
  ! Coeff3
  if (imod2==1) then
     coeff3(:,:)=efvm_cat0_coeffs2(:,:)
  else
     call reconstruct_efvm_mlt_coeffs(efvm_all_ca_coeffs2(imod1,:,:,:),nmlat,&
          LMAX,NMAX,ca_in,coeff3)
  end if

  ! Coeff4
  if (imod2==ncf) then
     coeff4(:,:)=coeff3(:,:)
  else
     call reconstruct_efvm_mlt_coeffs(efvm_all_ca_coeffs2(imod2,:,:,:),nmlat,&
          LMAX,NMAX,ca_in,coeff4)                                              
  end if

  coeff2(:,:) = coeff3(:,:) * wgt1 + coeff4(:,:) * wgt2

end subroutine calc_efvm_coeffs

!! ONE WAY TO CALCULATE DED1 AND DED2 WHEN THE CF IS LARGE
! -----------------------------------------------------------------------------
subroutine calc_extrapolated_dEd(coeff,LMAX,nmlat,nmlt,mlats,mlts,&
     er,disp,dEd_out)

  implicit none 

  integer, intent(in) :: LMAX, nmlat, nmlt
  real, intent(in) :: coeff(nmlat,2*LMAX+1), mlats(nmlat), mlts(nmlt)
  real, intent(in) ::  er, disp
  real, intent(out) :: dEd_out(nmlt,nmlat)

  real :: rs1(nmlat), ts1(nmlt), wgt1, wgt2, val
  integer :: imlt, imlat, ir1, ir2, kk
  real :: ts_1, rs_1, x1, y1, dx, dy, dx1, dy1, x2, y2, ts_2, rs_2
  real :: PARM(2*LMAX+1), fs(2*LMAX+1)

  real, parameter :: pi = 3.14159265359 

  ts1(:)=mlts(:)/12*pi
  rs1(:)=90-mlats(:)

  do imlt=1,nmlt
     do imlat=1,nmlat
  
        ts_1=ts1(imlt)
        rs_1=rs1(imlat)
        x1=rs_1*cos(ts_1)
        y1=rs_1*sin(ts_1)

        dx1=x1/er
        dy1=y1/er

        x2=dx1-disp
        y2=dy1
        
        rs_2=sqrt(x2**2+y2**2)
        ts_2=atan2(y2,x2)
        
        if (ts_2<0.) ts_2=ts_2+2*pi

        if (rs_2<=rs1(nmlat)) then
           ir1=nmlat
           ir2=nmlat
           wgt1=1.
           wgt2=0.
        else if (rs_2>rs1(1)) then
           ir1=0
           ir2=0
        else
           do kk=1,nmlat-1
              if ((rs1(kk)>=rs_2) .and. (rs1(kk+1)<rs_2)) then
                 ir1=kk
                 ir2=kk+1
                 wgt1=(rs_2-rs1(kk+1))/(rs1(kk)-rs1(kk+1))
                 wgt2=1-wgt1
                 exit 
              end if
           end do
        end if

        ! Reconstruct the Ed1
        if (ir1>0) then
           PARM(:)=coeff(ir1,:)*wgt1+coeff(ir2,:)*wgt2
           call construct_fourier_series(ts_2,LMAX,fs)
           val=dot_product(PARM,fs)
        else
           val=0.
        end if

        if (val<0.) val=0
        dEd_out(imlt,imlat) = val

     end do ! MLT
  end do ! MLAT

end subroutine calc_extrapolated_dEd
