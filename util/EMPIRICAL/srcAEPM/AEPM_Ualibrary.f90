! Higher-class functions of the AEPM
! Updated from the ModSAM.f90
! Created: Qingyu Zhu, 05/18/2020
!
! -----------------------------------------------------------------------------

!!! Reconstruct the MLT coeffs from the Ca coeffs if Ca_in is known 
! -----------------------------------------------------------------------------
subroutine reconst_aepm_mlt_coeffs(coeff_in,nchannel,nmlat,LMAX,NMAX,ca_in,&
     coeff_out)

  implicit none 

  integer, intent(in) :: nchannel,nmlat,LMAX,NMAX
  real, intent(in) :: coeff_in(nchannel,nmlat,2*LMAX+1,2*NMAX+1)
  real, intent(in) :: ca_in
  real, intent(out) :: coeff_out(nchannel,nmlat,2*LMAX+1)

  integer :: ichannel, imlat, iparm 
  real :: PARM(2*NMAX+1), fs(2*NMAX+1) 

  do ichannel=1,nchannel
     do imlat=1,nmlat
        do iparm=1,2*LMAX+1

           PARM(:)=coeff_in(ichannel,imlat,iparm,:)                         
           call construct_fourier_series(ca_in,NMAX,fs)                     
           coeff_out(ichannel,imlat,iparm)=dot_product(PARM,fs)

        end do
     end do
  end do

end subroutine reconst_aepm_mlt_coeffs

!!! Reconstruct the distribution of the diff_ef at different channels 
! -----------------------------------------------------------------------------
subroutine reconst_aepm_diff_ef(coeff_in,nchannel,nmlt,mlts,nmlat,LMAX,&
     diff_ef_out)

  implicit none 

  integer, intent(in) :: nchannel,nmlat,LMAX,nmlt
  real, intent(in) :: mlts(nmlt)
  real, intent(in) :: coeff_in(nchannel,nmlat,2*LMAX+1) 
  real, intent(out) :: diff_ef_out(nmlt,nmlat,nchannel)

  integer :: ichannel, imlat, imlt                                          
  real :: PARM(2*LMAX+1), fs(2*LMAX+1), phi, val

  real, parameter :: pi = 3.14159265359

  do ichannel=1,nchannel                                                
                                                                               
     do imlat=1,nmlat                                                  
        
        PARM=coeff_in(ichannel,imlat,:)                                    
        
        do imlt=1,nmlt                                                 
           
           phi=mlts(imlt)/12*pi                                        
           call construct_fourier_series(phi,LMAX,fs)                      
           
           val=dot_product(PARM,fs)                                        
           
           ! Avoid negative values                                         
           if (val<=0) val=1.0e-6                                         
           
           diff_ef_out(imlt,imlat,ichannel)=val                            
           
        end do
        
     end do
     
  end do

end subroutine reconst_aepm_diff_ef

!!! Obtain the slope and yint for a given Ca_in
! -----------------------------------------------------------------------------
subroutine reconst_aepm_slope_yint(coeff1_in,coeff2_in,nchannel,NMAX,ca_in,&
     slope,yint)

  use aepm_interface

  implicit none                                                                
                                                                               
  integer, intent(in) :: nchannel,NMAX
  real, intent(in) :: coeff1_in(nchannel,2*NMAX+1), &
       coeff2_in(nchannel,2*NMAX+1), ca_in
  real, intent(out) :: slope(nchannel), yint(nchannel)

  integer :: ichannel, iparm 
  real :: PARM(2*NMAX+1), fs(2*NMAX+1)                                         
                                                                               
  do ichannel=1,11                                                       
     
     ! Slope 
     PARM(:)=coeff1_in(ichannel,:)                         
     call construct_fourier_series(ca_in,NMAX,fs)                        
     slope(ichannel)=dot_product(PARM,fs)  
     
     ! Yint
     PARM(:)=coeff2_in(ichannel,:)   
     call construct_fourier_series(ca_in,NMAX,fs)                              
     yint(ichannel)=dot_product(PARM,fs)

  end do

  slope=slope*(1e-4)

  if (aepm_debug) then
     write(*,*) "================= Slope, yint ====================="
     write (*,*) slope
     write (*,*) yint
     write (*,*)
  end if

end subroutine reconst_aepm_slope_yint

!!! Correct the hemispheric integrated differential energy flux 
! -----------------------------------------------------------------------------
subroutine correct_hemi_int_diff_ef(cf_in,slope,yint,nchannel,int_diff_ef)

  implicit none 

  integer, intent(in) :: nchannel 
  real, intent(in) :: cf_in
  real, intent(in) :: slope(nchannel), yint(nchannel)
  real, intent(out) :: int_diff_ef(nchannel)

  integer :: ichannel 

  int_diff_ef=0.
  
  do ichannel=1,11
     int_diff_ef(ichannel)=slope(ichannel)*cf_in+yint(ichannel)
  end do

end subroutine correct_hemi_int_diff_ef

!!! Calculate the expansion 
! -----------------------------------------------------------------------------
subroutine calc_aepm_expansion(cf_in,cf_ref,ca_in,cf_inf,expansion_rate)

  implicit none 

  real, intent(in) :: cf_in, cf_ref, ca_in, cf_inf
  real, intent(out) :: expansion_rate 

  real :: beta_in, beta_ref
  real :: half_theta, sin_ht, ca_factor
  real, parameter :: pi = 3.14159265359

  ! Calculate betas 
  call calc_beta(cf_in,cf_inf,beta_in)
  call calc_beta(cf_ref,cf_inf,beta_ref)

  half_theta=(ca_in)/2.
  sin_ht=sin(half_theta)
  ca_factor=sin_ht**2

  expansion_rate=1+ca_factor*&
       ((6.65*(1e-4)*beta_in+10.45)/(6.65*(1e-4)*beta_ref+10.45)-1)

end subroutine calc_aepm_expansion

!!! Calculate diff_ef at a hemisphere 
! -----------------------------------------------------------------------------
subroutine calc_aepm_diff_ef(cf_in,ca_in,rec_diff_ef)

  use aepm_interface

  implicit none 

  real, intent(in) :: cf_in, ca_in
  real, intent(out) :: rec_diff_ef(aepm_nmlt,aepm_nmlat,aepm_nchannel)

  integer :: imod1, imod2 
  real :: wgt1, wgt2 

  integer :: istat, ichannel,imlat 
  real :: mlats(aepm_nmlat+1), diff_cf, scale, expansion_rate 
  real, dimension(aepm_nchannel) :: rate, int_diff_ef1, &
       int_diff_ef2, int_diff_ef3, int_diff_ef4
  
  real, allocatable, dimension(:,:,:) :: mlt_coeff1, mlt_coeff2, rec_diff_ef1
  real, allocatable, dimension(:,:,:) :: diff_ef1, diff_ef2

  expansion_rate=1. 
  rec_diff_ef=0.

  !! Determine the modes and weights 
  call determine_cf_models(cf_in,aepm_ncf,aepm_ref_cfs,imod1,imod2)
  call determine_weights(cf_in,aepm_cf0,aepm_ncf,aepm_ref_cfs,&
       imod1,imod2,wgt1,wgt2)

  if (aepm_debug) &
       write (*,*) "Models and Weights: ", imod1, imod2, wgt1, wgt2

  !! Reconstruct the two modes 
  ! Diff_ef_2 

  if (.not. allocated(diff_ef1)) then 
     allocate(diff_ef1(aepm_nmlt,aepm_nmlat,aepm_nchannel),stat=istat)
  end if

  if (aepm_debug) &
       write (*,*) "istat=",istat

  if (imod2==1) then
     
     call reconst_aepm_diff_ef(aepm_cat0_coeffs,aepm_nchannel,&
          aepm_nmlt,aepm_mlts,aepm_nmlat,aepm_LMAX,diff_ef1)

  else
     
     if (.not. allocated(mlt_coeff1)) then
        allocate(mlt_coeff1(aepm_nchannel,aepm_nmlat,2*aepm_LMAX+1),&
             stat=istat)
     end if

     call reconst_aepm_mlt_coeffs(aepm_all_ca_coeffs(imod1,:,:,:,:),&
          aepm_nchannel,aepm_nmlat,aepm_LMAX,aepm_NMAX,ca_in,&    
          mlt_coeff1)

     call reconst_aepm_diff_ef(mlt_coeff1,aepm_nchannel,&
          aepm_nmlt,aepm_mlts,aepm_nmlat,aepm_LMAX,diff_ef1)
     
     deallocate(mlt_coeff1)

  end if

  if (aepm_debug) write (*,*) "Done calculating diff_ef1"

  ! Diff_ef_2
  if (.not. allocated(diff_ef2)) then                                          
     allocate(diff_ef2(aepm_nmlt,aepm_nmlat,aepm_nchannel),stat=istat)         
  end if

  if (aepm_debug) &
       write (*,*) "istat=",istat

  if (imod1==aepm_ncf) then
     diff_ef2=diff_ef1
  else
     if (.not. allocated(mlt_coeff2)) then
        allocate(mlt_coeff2(aepm_nchannel,aepm_nmlat,2*aepm_LMAX+1),&
             stat=istat)
     end if

     call reconst_aepm_mlt_coeffs(aepm_all_ca_coeffs(imod2,:,:,:,:),&
          aepm_nchannel,aepm_nmlat,aepm_LMAX,aepm_NMAX,ca_in,&    
          mlt_coeff2)

     call reconst_aepm_diff_ef(mlt_coeff2,aepm_nchannel,&
          aepm_nmlt,aepm_mlts,aepm_nmlat,aepm_LMAX,diff_ef2)

     deallocate(mlt_coeff2)

  end if

  if (aepm_debug) write (*,*) "Done calculating diff_ef2"

  rec_diff_ef=diff_ef1*wgt1+diff_ef2*wgt2

  deallocate(diff_ef1)
  deallocate(diff_ef2)

  !! Calcuate the hemispheric integrated differential energy flux 
  mlats = (/(imlat, imlat = 1,aepm_nmlat+1)/)+49. 

  call calc_hemi_int_val_channel(aepm_nmlt,aepm_nmlat,aepm_nchannel,&          
           rec_diff_ef,mlats,1.0e-6,1.,int_diff_ef1)

  call reconst_aepm_slope_yint(aepm_slope_coeffs,aepm_yint_coeffs,&
       aepm_nchannel,aepm_NMAX,ca_in,aepm_slope,aepm_yint)

  call correct_hemi_int_diff_ef(cf_in,aepm_slope,aepm_yint, &
       aepm_nchannel,int_diff_ef2)

  if (.not. allocated(rec_diff_ef1)) then                                     
     allocate(rec_diff_ef1(aepm_nmlt,aepm_nmlat,aepm_nchannel),stat=istat) 
  end if 

  rec_diff_ef1=rec_diff_ef

  ! Scale the diff_ef to the estimated values for the first 11 channels  
  do ichannel=1,11

     if (int_diff_ef1(ichannel)<=0) then 
        scale=1.
     else
        scale=int_diff_ef2(ichannel)/int_diff_ef1(ichannel)
     end if

     if (aepm_debug) &
          write (*,*) int_diff_ef1(ichannel), int_diff_ef2(ichannel), scale
     if (aepm_debug) write (*,*)

     rec_diff_ef1(:,:,ichannel)=rec_diff_ef(:,:,ichannel) * scale 

  end do
   
  rec_diff_ef(:,:,:) = rec_diff_ef1(:,:,:)

  if (aepm_debug) write (*,*) "Done correcting >500 eV diff_ef"
  if (aepm_debug) write (*,*) rec_diff_ef(:,15,4)

  aepm_expansion_rate=expansion_rate

  !! For large CFs 
  ! ---------------------------------------------------------------------------
  if (imod1==aepm_ncf) then

     if (.not. allocated(diff_ef1)) then 
        allocate(diff_ef1(aepm_nmlt,aepm_nmlat,aepm_nchannel),stat=istat)
     end if

     
     if (.not. allocated(mlt_coeff1)) then
        allocate(mlt_coeff1(aepm_nchannel,aepm_nmlat,2*aepm_LMAX+1),&
             stat=istat)
     end if

     imod1=imod1-1
     call reconst_aepm_mlt_coeffs(aepm_all_ca_coeffs(imod1,:,:,:,:),&
          aepm_nchannel,aepm_nmlat,aepm_LMAX,aepm_NMAX,ca_in,&    
          mlt_coeff1)

     call reconst_aepm_diff_ef(mlt_coeff1,aepm_nchannel,&
          aepm_nmlt,aepm_mlts,aepm_nmlat,aepm_LMAX,diff_ef1)
     
     deallocate(mlt_coeff1)

     ! Extrapolations for <0.5 keV channels 
     call calc_hemi_int_val_channel(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
          diff_ef1,mlats,1.0e-6,1.,int_diff_ef1)

     call calc_hemi_int_val_channel(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
          rec_diff_ef,mlats,1.0e-6,1.,int_diff_ef2)

     rate=(int_diff_ef2-int_diff_ef1)/&
          (aepm_ref_cfs(aepm_ncf)-aepm_ref_cfs(aepm_ncf-1))

     diff_cf = cf_in - aepm_ref_cfs(aepm_ncf)

     if (diff_cf<0) diff_cf=0. 

     int_diff_ef3=int_diff_ef2+rate*diff_cf

     do ichannel=12,aepm_nchannel

        if (int_diff_ef2(ichannel)<=0) then 
           scale=1.
        else
           scale=int_diff_ef3(ichannel)/int_diff_ef2(ichannel)
        end if

        rec_diff_ef1(:,:,ichannel)=rec_diff_ef(:,:,ichannel) * scale 

     end do

     ! Expansion
     call calc_aepm_expansion(cf_in,aepm_ref_cfs(aepm_ncf),ca_in,&
          aepm_cf_inf,expansion_rate)
    
     aepm_expansion_rate=expansion_rate
     
     call calc_hemi_int_val_channel(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
          rec_diff_ef1,mlats,1.0e-6,1.,int_diff_ef3)

     call calc_hemi_int_val_channel(aepm_nmlt,aepm_nmlat,aepm_nchannel,&
          rec_diff_ef1,mlats,1.0e-6,expansion_rate,int_diff_ef4)

     if (aepm_debug) write(*,*) "int_diff_ef4", int_diff_ef4

     ! Scale back from int_diff_ef4 to int_diff_ef3
     do ichannel=1,aepm_nchannel

        if (int_diff_ef4(ichannel)<=0) then 
           scale=1.
        else
           scale=int_diff_ef3(ichannel)/int_diff_ef4(ichannel)
        end if

        rec_diff_ef(:,:,ichannel)=rec_diff_ef1(:,:,ichannel) * scale 

     end do
     
  end if

  ! Get the expanded grids 
  call calc_expanded_grid(aepm_nmlat,aepm_mlats,&
       expansion_rate,aepm_mlats1)

  if (aepm_debug) write(*,*) "Expansion_rate:", expansion_rate
  if (aepm_debug) write(*,*) aepm_mlats1

  deallocate(rec_diff_ef1)

end subroutine calc_aepm_diff_ef
  
