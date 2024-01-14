! MODULE THAT CAN USE FAC TO CALCULATE HIGH-LATITUDE ELECTRIC POTENTIAL
! AIMED TO HANDLE ASYMMETRIC FACS
! UPDATED FROM MODFAC.F90
!
! -----------------------------------------------------------------------------
module ModFAC

  use params_module, only: nmlon, nmlat_h, nmlat_T1, nlonlat
  use coef_module, only: coef_ns2, coef_new1, lhs_ns, rhs_ns, cpcp_nh, cpcp_sh

  implicit none 

  real :: pfrac_fac(nmlat_h)
  real, parameter :: ZigP_min=1.5

  contains

    ! STEP 1: Calculate potential associated with wind dynamo
    ! -------------------------------------------------------------------------
    subroutine calc_dynamopotential

      ! Two assumptions made in this routine:
      ! 1) Radial currents generated through wind dynamo do not contribute
      !    to the FAC measured by AMPERE (check with Art)     
      ! 2) Wind dynamo is symmetric                                         

      implicit none

      write (*,*) "Done calculate wind dynamo potential"

    end subroutine calc_dynamopotential


    ! Step 2: Correct FAC 
    ! -------------------------------------------------------------------------
    subroutine correct_FAC

      ! Currently, we use -sum(FAC*area)/sum(|FAC*area|) as the corrector
      ! In the future, Pedersen correction can be tested

      use fieldline_p_module,only:  fieldline_p,fline_p
      use fieldline_s_module,only:  fieldline_s1,fline_s1 

      real, dimension(nmlon,nmlat_h,2) :: fac_hl
      real :: sum,sumP,corr,sumn,sums,area_fac
      real :: sumn_af, sums_af, corr_n, corr_s

      integer :: i, j, isn

      fac_hl=0.                                                                
      sum=0.                                                                   
      sumn=0.                                                                  
      sums=0.

      sums_af=0.                                                               
      sumn_af=0.

      ! Calculate the correction factor
      do isn=1,2                                                               
         do i=1,nmlon                                                          
            do j=2,nmlat_h                                                     
                                                                            
               ! Jr needs to be mapped to 110 km
               ! Therefore a factor of (r2/r1)^3*(sinI2/sinI1) is needed
               ! 1 --> 110 km on the FAC, r1=6481.0 km 
               ! 2 --> 780 km on the FAC, r2=7151.0 km
               ! (r2/r1)^3 = 1.343
               ! (sinI2/sinI1) is typically <1.02, thus is ignored 
               area_fac = fline_p(i,j,isn)%M3(1)                               
               fac_hl(i,j,isn) = fline_p(i,j,isn)%fac_hl1 * area_fac * 1.343

               if (fline_s1(i,j,isn)%zigP<ZigP_min) fac_hl(i,j,isn)=0.

               if (fac_hl(i,j,isn) .ne. 0) then
                  if(isn.eq.1) sums= sums+fac_hl(i,j,isn)    
                  if(isn.eq.1) sums_af= sums_af+abs(fac_hl(i,j,isn))
                  if(isn.eq.2) sumn= sumn+fac_hl(i,j,isn)
                  if(isn.eq.2) sumn_af= sumn_af+abs(fac_hl(i,j,isn)) 
               end if

            end do
         end do
      end do

      corr_n = sumn/sumn_af 
      corr_s = sums/sums_af 
      write(*,*) "Corrections for NH and SH are: ", corr_n, corr_s

      sum=0.
      sumn=0.                                                                  
      sums=0.

      ! Correct FAC
      do isn=1,2 

         if (isn==1) corr=corr_s    
         if (isn==2) corr=corr_n 
         
         do i=1,nmlon                                                          
            do j=2,nmlat_h 

               if (fac_hl(i,j,isn).ne.0) &
                    fac_hl(i,j,isn) = fac_hl(i,j,isn)-&
                    abs(fac_hl(i,j,isn))*corr 

               if(isn.eq.1) sums= sums+fac_hl(i,j,isn) 
               if(isn.eq.2) sumn= sumn+fac_hl(i,j,isn) 
               
               area_fac = fline_p(i,j,isn)%M3(1)

               if (area_fac>0) then                                            
                  fline_p(i,j,isn)%fac_hl = fac_hl(i,j,isn)/area_fac           
               else                                                            
                  fline_p(i,j,isn)%fac_hl = 0. 
               end if
               
            end do
         end do
      end do
      
      write(6,*) "Total NH and SH FAC after correction:", sumn, sums  

    end subroutine correct_FAC

    ! Step 3: Set a lower boundary for potential transition 
    ! -------------------------------------------------------------------------
    subroutine set_fac_hemis

      ! Notes by Qingyu Zhu, 03/29/2021
      ! Previous experiments suggested that there are always unrealistic
      !    structures below the region where FAC is zero;
      ! This may because the transition is not well set;  
      ! Previous a fixed transition zone at 40-45 MLAT is set, which may not 
      !    be ideal since FAC typically terminated at > 55 MLAT;
      ! Here, we will try to set a dynamic boundary based on the FAC
      !    distribution: perhaps we could use either lowest MLAT at all
      !    longitudes or lowest MLAT at different longitudes;              
      ! This time, we use the former boundary;
      ! May need to remove suspicious low-latitude FAC structure; 
    
      use params_module, only: pi
      use fieldline_p_module,only: fline_p  

      integer :: isn, i, j 
      real, parameter :: crit1 = 45., crit2 = 40.  
      real :: mlatin                           
         
      ! Hemisphere 
      do isn=1,1

         ! Set up pfrac_fac in each hemisphere 
         do j=1,nmlat_h                                                    

            mlatin=fline_p(1,j,isn)%mlat_m/pi*180.                          
            mlatin=abs(mlatin)

            if (mlatin<=crit2) then                                            
               pfrac_fac(j) = 0.                          
            else if (mlatin>=crit1) then                                       
               pfrac_fac(j) = 1.                          
            else                                                               
               pfrac_fac(j) = (mlatin-crit2)/(crit1-crit2)                
            end if

         end do                                                                

      end do

    end subroutine set_fac_hemis

    ! Step 4: calculate potential associated with FAC in each hemisphere 
    ! -------------------------------------------------------------------------
    subroutine calc_hemisphere_pot(iHemisphere) 

      use params_module, only: dtr 
      use fieldline_p_module,only: fline_p
      use fieldline_s_module,only: fline_s1
      use readhlpoten_module, only: poten_hl

      implicit none                                                            
      
      integer, intent(in) :: iHemisphere                                       
      integer :: isn, i, j, k, ic, istat, jj     
      real :: area_fac, dlatm
      real :: sum1, sum_af, corr

      real, allocatable :: pot_out(:,:)       

      ! Recalculate FAC
      integer :: im, ip
      real, allocatable :: fac_out(:,:)
      
      if (.not. allocated(coef_new1)) &
           allocate(coef_new1(nmlon,nmlat_h,10),stat=istat)
      
      coef_new1=0.

      isn = iHemisphere
      
      do i=1,nmlon
         do j=1,nmlat_h  

            area_fac = fline_p(i,j,isn)%M3(1)
             
            do ic=1,9      
               coef_new1(i,j,ic) = coef_ns2(i,j,iHemisphere,ic)  
            end do
            
            ! qingyu, 05/26/2021
            ! It seems that without the neutral wind dynamo on the RHS
            ! the cross-track ion drift as well as the Joule heating 
            ! can be significantly reduced (!)
            ! Therefore the neutral wind dynamo is added on the RHS
            ! Since it is assumed that the wind dynamo is symmetric 
            ! below 50 MLAT so it is reasonable to use asymmetric pattern
            ! At high latitudes

            coef_new1(i,j,10) = fline_p(i,j,iHemisphere)%fac_hl * area_fac &
                 + coef_ns2(i,j,iHemisphere,10)
            
            ! Correction 
            do ic=1,8  
               coef_new1(i,j,ic) = coef_new1(i,j,ic) * pfrac_fac(j)
            end do

            dlatm=0.032725 ! Use TIEGCM values

            coef_new1(i,j,9) = coef_new1(i,j,9) * pfrac_fac(j) + & 
                 (1-pfrac_fac(j)) * coef_new1(i,j,9) * & 
                 (dlatm/(10.*dtr))**2 

            coef_new1(i,j,10) = coef_new1(i,j,10) * pfrac_fac(j)

         end do
      end do
      
      if (.not. allocated(pot_out)) &
           allocate(pot_out(nmlon,nmlat_h),stat=istat) 
      
      pot_out=0
      call pot_solver(coef_new1,pot_out) 

      ! -----------------------------------------------------------------------
      ! Recalculate FAC 
      if (.not. allocated(fac_out)) & 
           allocate(fac_out(nmlon,nmlat_h),stat=istat)

      fac_out=0.
      
      isn=iHemisphere 

      do i=1,nmlon
         if (i==1) then
            im=nmlon
         else
            im=i-1
         end if
         
         if (i==nmlon) then 
            ip=1
         else
            ip=i+1
         end if


         do j=2,nmlat_h-1


            fac_out(i,j)=pot_out(ip,j)*coef_new1(i,j,1) + &
                 pot_out(ip,j+1)*coef_new1(i,j,2) + &
                 pot_out(i ,j+1)*coef_new1(i,j,3) + & 
                 pot_out(im,j+1)*coef_new1(i,j,4) + &
                 pot_out(im,j  )*coef_new1(i,j,5) + &
                 pot_out(im,j-1)*coef_new1(i,j,6) + &
                 pot_out(i ,j-1)*coef_new1(i,j,7) + & 
                 pot_out(ip,j-1)*coef_new1(i,j,8) + &
                 pot_out(i ,j  )*coef_new1(i,j,9) - &
                 coef_ns2(i,j,iHemisphere,10)


            area_fac = fline_p(i,j,isn)%M3(1) 

            if (area_fac .ne. 0) &
                 fline_p(i,j,isn)%fac_hl2=fac_out(i,j)/area_fac


         end do
      end do

      ! -----------------------------------------------------------------------

      ! Remove the boundary value (for neutral dynamo)
      pot_out=pot_out-pot_out(50,nmlat_h)

      isn=iHemisphere
      
      do i=1,nmlon                                                           
         do j=1,nmlat_h                     

            if (isn.eq.1) then                                                 
               jj=j                                                            
            else                                                               
               jj=nmlat_T1 - j + 1                                             
            end if

            poten_hl(i,jj) = pot_out(i,j)
            
            fline_p(i,j,isn)%pot = pot_out(i,j)
            fline_p(i,j,isn)%pot1 = pot_out(i,j)

         end do
      end do

      ! Calculate the CPCP in each hemisphere
      ! qingyu, 05/23/2021
      if (isn .eq. 1) cpcp_sh=maxval(pot_out)-minval(pot_out)                  
      if (isn .eq. 2) cpcp_nh=maxval(pot_out)-minval(pot_out)            

      if (allocated(pot_out)) deallocate(pot_out) 
      if (allocated(coef_new1)) deallocate(coef_new1)     
      
    end subroutine calc_hemisphere_pot
     
    ! Compact routine 
    ! -------------------------------------------------------------------------
    subroutine get_asym_pot

      implicit none 

      call correct_fac
      call set_fac_hemis  
      call calc_hemisphere_pot(2)
      call calc_hemisphere_pot(1)

      ! Do not deallocate coef_ns2 
      !if (allocated(coef_ns2)) deallocate(coef_ns2)

    end subroutine get_asym_pot

    ! Solver (qingyu, 11/26/2020)
    ! Adapted from subroutine const_rhs
    ! -------------------------------------------------------------------------
    subroutine pot_solver(coef_in,pot_out)
      
      use fieldline_p_module, only: fline_p
      
      implicit none
      
      real, intent(in) :: coef_in(nmlon,nmlat_h,10)
      real, intent(out) :: pot_out(nmlon,nmlat_h)
      
      integer :: i,j,im,ip,ijm,it,ij,jt,istat,isn
      
      integer:: info,ifail
      integer,parameter :: nrhmax = 1
      integer ::  ipiv(nlonlat)
      character, parameter :: trans='N'
      
      real(8) :: colsum,anorm ,rcond,work(4*nlonlat) 
      real(8) :: dlange 
      integer :: info_c,iwork(nlonlat)
      character, parameter :: norm='1' 
      
      if (.not. allocated(lhs_ns)) then
         allocate(lhs_ns(nlonlat,nlonlat),STAT=istat)
         allocate(rhs_ns(nlonlat),STAT=istat)
      end if
      
      lhs_ns=0.
      rhs_ns=0.
      ij=0
      isn=1
      pot_out=0.
      
      do i=1,nmlon
         
         if (i==1) then
            im=nmlon
         else
            im=i-1
         end if
     
         if (i==nmlon) then
            ip=1
         else
            ip=i+1
         end if
         
         j=1
         ij = (i-1)*nmlat_h+j
         rhs_ns(ij) = coef_in(i,j,10)
         
         it = (ip-1)*nmlat_h+j
         lhs_ns(ij,it)  = coef_in(i,j,1)
         
         it = (ip-1)*nmlat_h+j+1 
         lhs_ns(ij,it)  = coef_in(i,j,2)
         
         it = (i-1)*nmlat_h+j+1
         lhs_ns(ij,it)  = coef_in(i,j,3)
         
         it = (im-1)*nmlat_h+j+1
         lhs_ns(ij,it)  = coef_in(i,j,4)
         
         it = (im-1)*nmlat_h+j
         lhs_ns(ij,it)  = coef_in(i,j,5)
     
         it = (i-1)*nmlat_h+j 
         lhs_ns(ij,it)  = coef_in(i,j,9)
         
         do j=2,nmlat_h-1
            
            ij = (i-1)*nmlat_h+j                                               
            rhs_ns(ij) = coef_in(i,j,10)                            
            
            it = (ip-1)*nmlat_h+j                               
            lhs_ns(ij,it)  = coef_in(i,j,1)             
            
            it = (ip-1)*nmlat_h+j+1                                            
            lhs_ns(ij,it)  = coef_in(i,j,2)                
            
            it = (i-1)*nmlat_h+j+1                                             
            lhs_ns(ij,it)  = coef_in(i,j,3)
            
            it = (im-1)*nmlat_h+j+1             
            lhs_ns(ij,it)  = coef_in(i,j,4)     
            
            it = (im-1)*nmlat_h+j                                              
            lhs_ns(ij,it)  = coef_in(i,j,5)                           
            
            it = (im-1)*nmlat_h+j-1                              
            lhs_ns(ij,it)  = coef_in(i,j,6)                          
            
            it = (i-1)*nmlat_h+j-1                                
            lhs_ns(ij,it)  = coef_in(i,j,7)                   
            
            it = (ip-1)*nmlat_h+j-1                         
            lhs_ns(ij,it)  = coef_in(i,j,8)                      
            
            it = (i-1)*nmlat_h+j                                            
            lhs_ns(ij,it)  = coef_in(i,j,9)
         end do
         
         j=nmlat_h                                                       
         ij = (i-1)*nmlat_h+j                                      
         rhs_ns(ij) = coef_in(i,j,10)                                
         
         it = (ip-1)*nmlat_h+j                                       
         lhs_ns(ij,it)  = coef_in(i,j,1)                                
         
         it = (im-1)*nmlat_h+j                                       
         lhs_ns(ij,it)  = coef_in(i,j,5)                              
         
         it = (im-1)*nmlat_h+j-1                                      
         lhs_ns(ij,it)  = coef_in(i,j,6)                        
     
         it = (i-1)*nmlat_h+j-1                                      
         lhs_ns(ij,it)  = coef_in(i,j,7)                       
         
         it = (ip-1)*nmlat_h+j-1                                  
         lhs_ns(ij,it)  = coef_in(i,j,8)                          
         
         it = (i-1)*nmlat_h+j                    
         lhs_ns(ij,it)  = coef_in(i,j,9)
         
      end do
      
      ! Solve the LHS X = RHS
      anorm = 0.d0
  
      do i=1,nlonlat
         colsum=0.d0
         
         ! Compute norm-1 of A
         do j=1,nlonlat
            colsum = colsum + abs(lhs_ns(j,i))
         end do
         
         anorm = max(anorm,colsum)
         
      end do
      write(6,*) 'anorm ',anorm
      
      write (6,*) 'factorize A'
      call DGETRF(nlonlat,nlonlat,lhs_ns,nlonlat,ipiv,info)
      write(6,*) 'info ' , info
      
      if (info .ne. 0) then
         write(*,*) "Singular matrix"
         write(*,*) "Error in LHS and RHS constructions, check !!!"
      else
         CALL DGECON(norm,nlonlat,lhs_ns,nlonlat,anorm,rcond,&
              work,iwork,info_c)
         write(6,*) 'info_c ' , info_c 
         write(6,*) 'condit=', 1/rcond
         
         write (6,*) 'solve'
         call DGETRS(trans,nlonlat,nrhmax,lhs_ns,nlonlat,ipiv,&
              rhs_ns,nlonlat,info)
         write (6,*) 'after solve',info 
         
         it=0
         do i=1,nmlon
            do j=1,nmlat_h
               
               it=it+1
               pot_out(i,j)=rhs_ns(it)
               
            end do
         end do
     
      end if
  
    end subroutine pot_solver

end module ModFAC
