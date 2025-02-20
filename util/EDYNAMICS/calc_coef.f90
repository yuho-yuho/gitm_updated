     module coef_module
     
     use params_module, only: nmlon,nmlat_h,nhgt_fix,nlonlat,nmlatS2_h,nmlat_T1
     
     implicit none
     
!     real :: coef(nmlon,nmlat_h,nhgt_fix,9,2) ! for each P-point each hemisphere and height
!     real :: coef_ns(nmlon,nmlat_h,10)        ! lhs+rhs: for each P-point 
!     real(8) :: rhs_ns(nlonlat)                  ! forcing: for each P-point 
!     real(8) :: lhs_ns(nlonlat,nlonlat)          ! lhs: for each P-point 
     
     real, allocatable :: coef(:,:,:,:,:)       ! for each P-point each hemisphere and height


     real, allocatable :: coef_ns(:,:,:)        ! lhs+rhs: for each P-point 
     real, allocatable :: coef_ns2(:,:,:,:)     ! lhs+rhs: for each P-point 

     ! qingyu, 11/25/2020 
     ! used to document NH and SH coefficients for asymmetric FACs
     real, allocatable :: coef_ns3(:,:,:,:), coef_ns4(:,:,:,:)

     ! qingyu, 11/26/2020
     real, allocatable :: coef_new1(:,:,:), coef_new2(:,:,:)

     ! qingyu, 01/27/2021
     real, allocatable :: sym_pot(:,:)

     ! qingyu, 01/22/2021
     real :: pfrac(nmlat_h) ! same for both hemisphere 

     real(8), allocatable :: rhs_ns(:)          ! forcing: for each P-point 
     real(8), allocatable :: lhs_ns(:,:)        ! lhs: for each P-point 
     
     ! qingyu, 05/23/2021  
     real :: cpcp_nh, cpcp_sh

     contains

!-----------------------------------------------------------------------------      
     subroutine calc_coef
!        
     use fieldline_s_module,only: fieldline_s1,fline_s1, &
          fieldline_s2,fline_s2
     use fieldline_p_module,only:  fieldline_p,fline_p
	  
     implicit none
!     
! coef ordering from tiegcm
! ^   equatorward
! coef(4) (i-1,j+1)      coef(3) (i,j+1)    coef(2) (i+1,j+1)
! coef(5) (i-1,j)        coef(9) (i,j)	    coef(1) (i+1,j)
! coef(6) (i-1,j-1)      coef(7) (i,j-1)    coef(8) (i+1,j-1)
! v   poleward
!
! relationship between P,S1, and S2 point for the same index (i,j)
!  P(i,j) then is really S1(i+0.5,j) and S2(i,j+0.5) with j increasing equatorward
!  coefficient is calculated at P points
!   

      integer :: isn,i,j,k,im,nmax,status,ic
      real :: N2p_p,N2h_p
      
! allocate arrays
      allocate(coef(nmlon,nmlat_h,nhgt_fix,9,2),STAT=status) ! for each P-point each hemisphere and height
      if(status /= 0 ) write(6,*) 'alloc coef failed'
      allocate(coef_ns(nmlon,nmlat_h,10),STAT=status)	 ! lhs+rhs: for each P-point 
      if(status /= 0 ) write(6,*) 'alloc coef_ns failed'
      allocate(rhs_ns(nlonlat),STAT=status) 		 ! forcing: for each P-point 
      if(status /= 0 ) write(6,*) 'alloc rhs_ns failed'
      allocate(lhs_ns(nlonlat,nlonlat),STAT=status) 	 ! lhs: for each P-point 
      if(status /= 0 ) write(6,*) 'alloc lhs_ns failed'
      allocate(coef_ns2(nmlon,nmlat_h,2,10),STAT=status)  ! lhs+rhs: for each P-point both hemispheres 
      if(status /= 0 ) write(6,*) 'alloc coef_ns failed'
      

       do isn = 1,2 ! loop over both hemisphere
         do i=1,nmlon ! loop over all longitudes
	   if(i == 1) then
	     im = nmlon
	   else
	     im = i-1
	   endif
	 
           do j=2,nmlat_h ! loop over all latitudes in one hemisphere NOT THE POLE
	      nmax = fline_p(i,j,isn)%npts ! maximum of points on fieldline
	      
              do k=1,nmax  ! longitudinal wrap around
		if((nmlat_h-j+1) == k) then ! top volume at equator
		  N2h_p = 0.
		  N2p_p = 0.
		else
		  N2h_p = fline_s2(i,j,isn)%N2h(k)  ! i,j+0.5
		  N2p_p = fline_s2(i,j,isn)%N2p(k)  ! i,j+0.5
		endif 
	        coef(i,j,k,2,isn) = -fline_s1(i,j,isn)%N1h(k) +  N2h_p 
	        coef(i,j,k,3,isn) = fline_s1(im,j,isn)%N1h(k) - &
	         fline_s1(i,j,isn)%N1h(k)+ N2p_p
	        coef(i,j,k,4,isn) = fline_s1(im,j,isn)%N1h(k) - N2h_p	
!		 
	        coef(i,j,k,5,isn) = fline_s1(im,j,isn)%N1p(k) + &
	           fline_s2(i,j-1,isn)%N2h(k)-N2h_p
	        coef(i,j,k,6,isn) = -fline_s1(im,j,isn)%N1h(k) + &
	           fline_s2(i,j-1,isn)%N2h(k)
	        coef(i,j,k,7,isn) = -fline_s1(im,j,isn)%N1h(k) + &
	           fline_s1(i,j,isn)%N1h(k)+fline_s2(i,j-1,isn)%N2p(k)
	        coef(i,j,k,8,isn) = fline_s1(i,j,isn)%N1h(k) - &
	           fline_s2(i,j-1,isn)%N2h(k)
	        coef(i,j,k,1,isn) = fline_s1(i,j,isn)%N1p(k) - &
	           fline_s2(i,j-1,isn)%N2h(k)+N2h_p
	        coef(i,j,k,9,isn) = -fline_s1(im,j,isn)%N1p(k) - &
	           fline_s1(i,j,isn)%N1p(k)-fline_s2(i,j-1,isn)%N2p(k) -&
	           N2p_p  
              end do  ! end height loop
           end do  ! end lat/fieldline loop
!    	  	     
         end do  ! end longitude loop
       end do ! end hemisphere loop
! 
! add the coefficients in height to get coefficients for each hemisphere 
! needed for calculating high latitude FAC      
! initialize     
     coef_ns2 = 0.
 !   
     do i=1,nmlon   ! loop over all longitudes
       do isn = 1,2   ! hemisphere loop
        j=1               ! polar values are set
        coef_ns2(i,j,isn,9)   = 0.5  ! later added together hemispheres to get one
        coef_ns2(i,j,isn,1:8) = 0.
        coef_ns2(i,j,isn,10)  = 0.  ! set potential at pole
        do j=2,nmlat_h ! loop over all latitudes in one hemisphere NOT THE POLE
     	  nmax = fline_p(i,j,1)%npts ! maximum of points on fieldline     	  
     	  do k=1,nmax  ! height loop

           !if (isnan(coef(i,j,k,1,isn))) &
           !     write(6,*) "NaN exisits in coef for point #",&
           !     i,j,k,"Ignore ..."

	    do ic = 1,9 ! 9-point stencil
                 
              if (.not. isnan(coef(i,j,k,ic,isn))) then          
                 coef_ns2(i,j,isn,ic) = coef_ns2(i,j,isn,ic)+ &
                      coef(i,j,k,ic,isn)
              end if

	    end do

            if (.not. isnan(fline_p(i,j,isn)%S(k))) then
               coef_ns2(i,j,isn,10) = coef_ns2(i,j,isn,10)+ &
                    fline_p(i,j,isn)%S(k)
            end if
     	  end do  ! end height loop
!	  


        !if ((i.eq.2) .and. (j.eq.30)) write(*,*) "Before add_coef_ns", &
        !     i,j,isn,coef_ns2(i,j,isn,:)


!        if (isnan(coef_ns2(2,30,isn,1))) then
!           write(6,*) "NaN exists in coef_ns2"
!           do k=1,nmax
!              do ic=1,9
!                 write(6,*) k, ic, coef(2,30,k,ic,isn)
!              end do
!              ic=10
!              write(6,*) k, ic, fline_p(2,30,isn)%S(k)
!           end do
!        end if

        end do  ! end lat/fieldline loop
       end do  ! end hemisphere loop
!     	 
     end do  ! end longitude loop
!     
     end subroutine calc_coef

     ! qingyu, 01/28/2021
     ! turn off wind dynamo before adding FACs or prescribing high-latitude 
     ! electric potential
     ! Only for |MLAT|<50
     subroutine close_winddynamo

       use fieldline_p_module, only: fieldline_p,fline_p
       use params_module, only: nmlon, nmlat_h, pi

       implicit none 

       integer :: i,isn,j,idx
       real :: mlatin

       !coef_ns2(:,:,:,10)=0.

       do j=1,nmlat_h
          
          mlatin=fline_p(1,j,1)%mlat_m/pi*180.
          
          if (abs(mlatin)<50) then
             idx=j
             exit
          end if

       end do

       coef_ns2(:,idx:nmlat_h,:,10)=0.

     end subroutine close_winddynamo

     ! qingyu, 01/23/21
     ! Add an alternate methodology to construct the electrodynamo solver
     ! Try to include the penetration electric field
     ! Adapted from the TIEGCM procedure
     ! ------------------------------------------------------------------------
     subroutine calc_coef1

       implicit none 
       
       call set_transition 
       call correct_coef_ns2

     end subroutine calc_coef1

     ! Set up the transition zone
     ! ------------------------------------------------------------------------
     subroutine set_transition                                                
                                                                          
       use params_module, only: pi  
       use fieldline_p_module,only: fline_p 
       use readhlpoten_module, only: poten_hl

       implicit none                                                       

       integer :: isn, i, j                                             
       real :: crit1, crit2                                               
       real :: mlatin, cpcp_in                                                
                                                                              
       crit1=55.                                               
       crit2=50.

       ! Use a dynamic transition zone
       ! Currently use the NH CPCP to determine the poleward boundary
       ! of the transition zone 
       ! qingyu, 05/25/2021

       cpcp_nh=maxval(poten_hl(:,nmlat_h:nmlat_T1))-&
            minval(poten_hl(:,nmlat_h:nmlat_T1)) 

       cpcp_sh=maxval(poten_hl(:,1:nmlat_h))-&
            minval(poten_hl(:,1:nmlat_h))        

       cpcp_in = cpcp_nh/1000.

       crit1= 90-(-3.8 + 8.48 * cpcp_in**(0.1875) + 5) ! Heelis formula

       if (crit1>75) crit1=75.
       if (crit1<65) crit1=65.
       crit2=crit1-15. 

       do j=1,nmlat_h

          mlatin=fline_p(1,j,1)%mlat_m/pi*180. 
          mlatin=abs(mlatin) 
          
          if (mlatin<=crit2) then                                          
             pfrac(j) = 1.                                                    
          else if (mlatin>=crit1) then                                       
             pfrac(j) = 0.                                                 
          else                                                             
             pfrac(j) = 1.-(mlatin-crit2)/(crit1-crit2)                      
          end if

       end do

     end subroutine set_transition

     ! Correct the coef_ns2
     ! ------------------------------------------------------------------------
     subroutine correct_coef_ns2 

       use params_module, only: dtr, pi
       use fieldline_p_module,only: fline_p 
       use readhlpoten_module, only: poten_hl 

       implicit none 

       integer :: isn, i, j, jj, ic
       real :: dlatm
       integer :: status

      ! qingyu, 01/27/2021
      if (.not. allocated(sym_pot)) then
         allocate(sym_pot(nmlon,nmlat_h),STAT=status)
         if(status/=0 ) write(6,*) 'alloc sym_pot failed'
      end if     
       
       
       do isn=1,2

          sym_pot=0.
          
          do i=1,nmlon                                                        
             do j=1,nmlat_h                                                   
          

                if (isn .eq. 1) then                                        
                   jj=j    
                   !jj = nmlat_T1 - j + 1 
                else                                                      
                   jj = nmlat_T1 - j + 1                            
                end if
                                                                    
                if (j==1) then                                              
                   dlatm=abs(fline_p(i,j,isn)%mlat_m-fline_p(i,j+1,isn)%mlat_m)
                else                                                           
                   dlatm=abs(fline_p(i,j-1,isn)%mlat_m-fline_p(i,j,isn)%mlat_m)
                end if

                dlatm=0.032725 ! Use TIEGCM values 

                !write(*,*) i,j,isn, dlatm

                ! 01/27/2021
                ! Symmetric potential 
                ! (sigma_R(N)*phi(N)+sigma_R(S)*phi(S))/(sigma_R(N)+sigma_R(S))
                ! TIEGCM assumes phi(N)==phi(S)

                ! Get the symmetric potential before corrections 
                !if (sum(coef_ns2(i,j,:,9)) .ne. 0) then
                !   sym_pot(i,j) = (coef_ns2(i,j,1,9)*poten_hl(i,j)+&
                !        coef_ns2(i,j,2,9)*poten_hl(i,nmlat_T1-j+1))/&
                !        sum(coef_ns2(i,j,:,9))
                !end if

                ! Directly use NH potential as SYM pot 
                sym_pot(i,j) = poten_hl(i,nmlat_T1-j+1)

                ! Coefficient 10  
                ! Penetration electric fields do exist but need validations
                ! Therefore the wind dynamo term is kept
                coef_ns2(i,j,isn,10) = pfrac(j)*coef_ns2(i,j,isn,10)+&
                     (1-pfrac(j))*coef_ns2(i,j,isn,9)*(dlatm/(10.*dtr))**2*&
                     sym_pot(i,j)

                ! Coefficient 9
                coef_ns2(i,j,isn,9)=pfrac(j)*coef_ns2(i,j,isn,9) + &
                     (1-pfrac(j))*coef_ns2(i,j,isn,9)*(dlatm/(10.*dtr))**2


                ! Coefficient 1-8
                do ic=1,8                                           
                   coef_ns2(i,j,isn,ic) = coef_ns2(i,j,isn,ic) * pfrac(j)    
                end do

             end do
          end do
       end do              

     end subroutine correct_coef_ns2
     !-------------------------------------------------------------------------
     subroutine calc_FAC
!    
     use readhlpoten_module, only: poten_hl  ! P points
     use fieldline_p_module,only:  fieldline_p,fline_p
     use fieldline_s_module,only:  fieldline_s1,fline_s1
!
     implicit none
!
     real, dimension(nmlon,nmlat_h,2) :: fac_hl 
     real :: sum,sumP,corr,sumn,sums
!
     integer :: isn,i,j,jj,icof,im,ip
!     
     real :: area_fac, sumn_af, sums_af, corr_n, corr_s

     fac_hl = 0.
     sum    = 0.
     sumn   = 0.
     sums   = 0.
     sumP   = 0.

     sumn_af=0.
     sums_af=0.

!     coef_ns2(:,:,:,1:8) = 0. 
     do isn = 1,2 ! loop over both hemisphere
       do i=1,nmlon ! loop over all longitudes
         if(i == 1) then ! wrap around in longitude
           im = nmlon
         else
           im = i-1	  
         endif
         if(i == nmlon) then
           ip = 1
         else
           ip = i+1
         endif
         do j=2,nmlat_h ! loop over all latitudes in one hemisphere not the pole (potential set later)

           ! Use new correction method
           area_fac = fline_p(i,j,isn)%M3(1)

           if(isn.eq.1) then  ! jj is latitude index from pole to pole
	    jj = j
	    fac_hl(i,j,isn) = poten_hl(ip,jj)*coef_ns2(i,j,isn,1)+ &
	     poten_hl(ip,jj+1)*coef_ns2(i,j,isn,2)+ &
	     poten_hl(i ,jj+1)*coef_ns2(i,j,isn,3)+ &
	     poten_hl(im,jj+1)*coef_ns2(i,j,isn,4)+ &
	     poten_hl(im,jj  )*coef_ns2(i,j,isn,5)+ &
	     poten_hl(im,jj-1)*coef_ns2(i,j,isn,6)+ &
	     poten_hl(i ,jj-1)*coef_ns2(i,j,isn,7)+ &
	     poten_hl(ip,jj-1)*coef_ns2(i,j,isn,8)+ &
	     poten_hl(i  ,jj )*coef_ns2(i,j,isn,9)
            
            if (fac_hl(i,j,isn) .ne. 0) then
               sums= sums+fac_hl(i,j,isn)
               sums_af= sums_af+abs(fac_hl(i,j,isn))
            end if
               
	   else
	    jj = nmlat_T1 - j + 1
	    fac_hl(i,j,isn) = poten_hl(ip,jj)*coef_ns2(i,j,isn,1)+ &
	     poten_hl(ip,jj-1)*coef_ns2(i,j,isn,2)+ &
	     poten_hl(i ,jj-1)*coef_ns2(i,j,isn,3)+ &
	     poten_hl(im,jj-1)*coef_ns2(i,j,isn,4)+ &
	     poten_hl(im,jj  )*coef_ns2(i,j,isn,5)+ &
	     poten_hl(im,jj+1)*coef_ns2(i,j,isn,6)+ &
	     poten_hl(i ,jj+1)*coef_ns2(i,j,isn,7)+ &
	     poten_hl(ip,jj+1)*coef_ns2(i,j,isn,8)+ &
	     poten_hl(i  ,jj )*coef_ns2(i,j,isn,9)

            if (fac_hl(i,j,isn) .ne. 0) then
               sumn= sumn+fac_hl(i,j,isn)
               sumn_af= sumn_af+abs(fac_hl(i,j,isn))
            end if

	   endif
            !sum  = sum  + fac_hl(i,j,isn) 
	    !sumP = sumP + fline_s1(i,j,isn)%zigP*abs(fac_hl(i,j,isn))
	    !sumP = sumP + abs(fac_hl(i,j,isn))
	    !if(isn.eq.1) sums= sums+fac_hl(i,j,isn) 
	    !if(isn.eq.2) sumn= sumn+fac_hl(i,j,isn) 
     	 end do  ! end lat/fieldline loop 
         !	       	   
       end do  ! end longitude loop
     end do ! end hemisphere loop

     corr_n = sumn/sumn_af                                           
     write(6,*) "------ NH Total FAC", sumn                                   
     write(6,*) "------ NH correction", corr_n                                
     ! SH correction                                                          
     corr_s = sums/sums_af                                                    
     write(6,*) "------ SH Total FAC", sums                                   
     write(6,*) "------ SH correction", corr_s
!     
! next two lines to make sure no high latitude forcing is included    
! Commented to let the FAC go
!     fac_hl= 0.
!     corr =  0.


!
!     write(6,*) 'sum fac_hl', sum, maxval(fac_hl),minval(fac_hl)
!     corr = sum/sumP
!     write(6,*) 'correction', corr
!     write(6,*) 'sumn', sumn, 'sums ',sums
!

    
     sum = 0.
     sumn = 0.
     sums = 0.

     do isn = 1,2 ! loop over both hemisphere

       if (isn==1) corr=corr_s 
       if (isn==2) corr=corr_n 

       do i=1,nmlon ! loop over all longitudes
         do j=2,nmlat_h ! loop over all latitudes in one hemisphere not the pole (potential set later)
	   
	   !fac_hl(i,j,isn) = fac_hl(i,j,isn)-fline_s1(i,j,isn)%zigP*abs(fac_hl(i,j,isn))*corr


           area_fac = fline_p(i,j,isn)%M3(1) 

           if (fac_hl(i,j,isn).ne.0) &                                      
                !fac_hl(i,j,isn) = fac_hl(i,j,isn)-area_fac*corr
                fac_hl(i,j,isn) = fac_hl(i,j,isn)-abs(fac_hl(i,j,isn))*corr
	   sum = sum   + fac_hl(i,j,isn) 
           if (isn==1) sums=sums+fac_hl(i,j,isn) 
           if (isn==2) sumn=sumn+fac_hl(i,j,isn)
	   ! put in coef-array  
	   coef_ns2(i,j,isn,10)	=   coef_ns2(i,j,isn,10) + fac_hl(i,j,isn)
	   !
!	   write(66,'(3(x,e15.8))') fline_p(i,j,isn)%mlon_qd(1), &
!           fline_p(i,j,isn)%mlat_qd(1), coef_ns2(i,j,isn,10)  
	   if(fline_p(i,j,isn)%M3(1).ne.0) then
!	     write(77,'(3(x,e15.8))') fline_p(i,j,isn)%mlon_qd(1), &
!           fline_p(i,j,isn)%mlat_qd(1), coef_ns2(i,j,isn,10)/ fline_p(i,j,isn)%M3(1) 
!
	   !fline_p(i,j,isn)%fac_hl = coef_ns2(i,j,isn,10)/ fline_p(i,j,isn)%M3(1) 
       
       fline_p(i,j,isn)%fac_hl = fac_hl(i,j,isn)/area_fac 

	   endif
     	 end do  ! end lat/fieldline loop 
         !	       	   
       end do  ! end longitude loop
     end do ! end hemisphere loop
!     
     !write(6,*) 'after corr. sum fac_hl', sum, sums, sumn
	 
     !coef_ns2(:,:,1,:) = 0. 
!     
     end subroutine calc_FAC
!----------------------------------------------------------------------------- 

     subroutine add_coef_ns
! adds the coefficients from both hemispheres together
!  output:  coef_ns    
     use fieldline_p_module,only:  fieldline_p,fline_p

     implicit none
!
     integer :: i,j,k,ic,nmax,status
!
! add hemispheres together and set equatorial boundary condition   
     do i=1,nmlon ! loop over all longitudes
       do j=1,nmlat_h ! loop over all latitudes in one hemisphere
	  do ic = 1,9 ! 9-point stencil
	    coef_ns(i,j,ic) = coef_ns2(i,j,1,ic)+coef_ns2(i,j,2,ic)
	  end do
	  coef_ns(i,j,10) = coef_ns2(i,j,1,10)+ coef_ns2(i,j,2,10)
!	  
       end do  ! end lat/fieldline loop
       j=nmlat_h                    ! set equatorial values (page 14 Art's notes)
       coef_ns(i,j,2) = 0.
       coef_ns(i,j,3) = 0.
       coef_ns(i,j,4) = 0.  
!     	 
     end do  ! end longitude loop

! check the coef_ns
     !write(*,*) "latitude",fline_p(2,31,1)%mlat_m/3.14*180.

     !write(*,*) "In add_coef_ns SH coef_ns2",coef_ns2(2,30,1,:)
     !write(*,*) "In add_coef_ns NH coef_ns2",coef_ns2(2,30,2,:)

!     
! check stability     
!     do i=1,nmlon ! loop over all longitudes
!       do j=1,nmlat_h ! loop over all latitudes in one hemisphere NOT THE POLE
!         if(coef_ns(i,j,1)*coef_ns(i,j,5)<0 ) then
!	   write(88,*) 'lon ',fline_p(i,j,2)%mlon_m, &
!	       fline_p(i,j,2)%mlat_m
!	 end if
!         if(coef_ns(i,j,3)*coef_ns(i,j,7)<0 ) then
!	   write(88,*) 'lat ',fline_p(i,j,2)%mlon_m, &
!	       fline_p(i,j,2)%mlat_m
!	 end if
!       end do  ! end lat/fieldline loop
!     end do  ! end longitude loop
     
     deallocate(coef,STAT=status)  
     if(status /= 0) write(6,*) 'deallocation of coef not succesful'  
     
     deallocate(coef_ns2,STAT=status)  
     if(status /= 0) write(6,*) 'deallocation of coef_ns2 not succesful'  
     
     end subroutine add_coef_ns
!----------------------------------------------------------------------------- 
     subroutine const_rhs
! construct matrix lhs &rhs  for solving with matlab 

     use fieldline_p_module, only: fieldline_p,fline_p
!     
     implicit none
!
     integer :: i,j,im,ip,ijm,it,ij,jt,status,isn
!  for solver
     integer:: info,ifail  
     integer,parameter :: nrhmax = 1        
     integer ::  ipiv(nlonlat)   
     character, parameter :: trans='N'
! condition number    
     real(8) :: colsum,anorm ,rcond,work(4*nlonlat)
     real(8) :: dlange
     integer :: info_c,iwork(nlonlat)
     character, parameter :: norm='1'
! am 1/2015 for testing      
     real :: pot_ns(nlonlat)  
     real :: z_ns(nlonlat)
!    
     lhs_ns = 0.
     rhs_ns = 0.
     ij = 0
     isn = 1

! am 1/2015 test by putting analytical potential into stencil and compare RHS
!     
     do i=1,nmlon       ! loop over all longitudes
!       
       if(i == 1) then ! wrap around in longitude
         im = nmlon
       else
         im = i-1	
       endif
       if(i == nmlon) then
         ip = 1
       else
         ip = i+1
       endif
       ! pole values
       j=1  ! rhs = 0. c9=1, and c(1:8) = 0 for the pole
       ij = (i-1)*nmlat_h+j
       rhs_ns(ij) = coef_ns(i,j,10)
       it = (ip-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns(i,j,1)
       it = (ip-1)*nmlat_h+j+1
       lhs_ns(ij,it)  = coef_ns(i,j,2)
       it = (i-1)*nmlat_h+j+1
       lhs_ns(ij,it)  = coef_ns(i,j,3)
       it = (im-1)*nmlat_h+j+1
       lhs_ns(ij,it)  = coef_ns(i,j,4)
       it = (im-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns(i,j,5)
       it = (i-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns(i,j,9)
       !
       pot_ns(ij) = fline_p(i,j,isn)%pot_test ! am 1/2015 test
       !
       do j=2,nmlat_h-1   ! pole to equator
           ij = (i-1)*nmlat_h+j
           rhs_ns(ij) = coef_ns(i,j,10)
           it = (ip-1)*nmlat_h+j
           lhs_ns(ij,it)  = coef_ns(i,j,1)
           it = (ip-1)*nmlat_h+j+1
           lhs_ns(ij,it)  = coef_ns(i,j,2)
           it = (i-1)*nmlat_h+j+1
           lhs_ns(ij,it)  = coef_ns(i,j,3)
           it = (im-1)*nmlat_h+j+1
           lhs_ns(ij,it)  = coef_ns(i,j,4)
           it = (im-1)*nmlat_h+j
           lhs_ns(ij,it)  = coef_ns(i,j,5)
           it = (im-1)*nmlat_h+j-1
           lhs_ns(ij,it)  = coef_ns(i,j,6)
           it = (i-1)*nmlat_h+j-1
           lhs_ns(ij,it)  = coef_ns(i,j,7)
           it = (ip-1)*nmlat_h+j-1
           lhs_ns(ij,it)  = coef_ns(i,j,8)
           it = (i-1)*nmlat_h+j
           lhs_ns(ij,it)  = coef_ns(i,j,9)
       !
           pot_ns(ij) = fline_p(i,j,isn)%pot_test  ! am 1/2015 test
       !
       end do  ! end lat/fieldline loop
       ! equator values
       j=nmlat_h
       ij = (i-1)*nmlat_h+j
       rhs_ns(ij) = coef_ns(i,j,10)
       it = (ip-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns(i,j,1)
       it = (im-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns(i,j,5)
       it = (im-1)*nmlat_h+j-1
       lhs_ns(ij,it)  = coef_ns(i,j,6)
       it = (i-1)*nmlat_h+j-1
       lhs_ns(ij,it)  = coef_ns(i,j,7)
       it = (ip-1)*nmlat_h+j-1
       lhs_ns(ij,it)  = coef_ns(i,j,8)
       it = (i-1)*nmlat_h+j
       lhs_ns(ij,it)  = coef_ns(i,j,9)
       !
       pot_ns(ij) = fline_p(i,j,isn)%pot_test  ! am 1/2015 test
       !
     end do  ! end longitude loop
!     
!   
     deallocate(coef_ns,STAT=status)  
     if(status /= 0) write(6,*) 'deallocation of coef_ns not succesful' 
!
! put pot_ns into the lhs lhs_ns -> lhs_ns*pot_ns= rhs_test
!   
       z_ns= 0. 
        z_ns = matmul( lhs_ns,pot_ns)
!       it=0
!       do i=1,nmlon
!         do j=1,nmlat_h
!          it = it +1
!          write(55,'(2(i4,x),3(x,e15.8))') i,j,z_ns(it),rhs_ns(it),fline_p(i,j,isn)%pot_test
!         enddo
!       enddo
!   
!   
!     do i=1,nlonlat	 ! loop over all longitudes
!       do j=1,nlonlat  ! pole to equator
!	 write(88,'(2(i4,x),1(x,e15.8))') i,j,lhs_ns(i,j)
!       end do  ! end lat/fieldline loop
!       write(99,'(1(i4,x),1(x,e15.8))') i,rhs_ns(i)
!     end do  ! end longitude loop

! solve the matrix LHS X = RHS
!
! calculate norm
      anorm = 0.d0
      do i=1,nlonlat
        colsum = 0.d0
        do j=1,nlonlat
      ! compute norm-1 of A-> max_j SUM_i abs(a_ij) 
           colsum = colsum + abs(lhs_ns(j,i))
        end do 
	anorm = max(anorm,colsum)
      end do 
      write(6,*) 'anorm ',anorm
!
!       Factorize A
!
   
      write (6,*) 'factorize A'
      call DGETRF(nlonlat,nlonlat,lhs_ns,nlonlat,ipiv,info)
      write(6,*) 'info ' , info
!
      if (info.EQ.0) then
!      
! calculate condition number      
        CALL DGECON( norm,nlonlat,lhs_ns,nlonlat,anorm,rcond,work,iwork,info_c)
        write(6,*) 'info_c ' , info_c
	write(6,*) 'condit=', 1/rcond
!
!        Compute solution
!
! am 1/2015 test
!         rhs_ns = z_ns 
!       it=0
!       do i=1,nmlon
!         do j=1,nmlat_h
!          it = it +1
!          if(j.eq.97.or.j.eq.98) rhs_ns(it) = z_ns(it) 
!         enddo
!       enddo
!       it=0
!       do i=1,nmlon
!         do j=1,nmlat_h
!          it = it +1
!          write(55,'(2(i4,x),3(x,e15.8))') i,j,z_ns(it),rhs_ns(it),fline_p(i,j,isn)%pot_test
!         enddo
!       enddo
!	 
         write (6,*) 'solve'
         call DGETRS(trans,nlonlat,nrhmax,lhs_ns,nlonlat,ipiv,rhs_ns,nlonlat,info)
         write (6,*) 'after solve',info
!
!        Print solution
!
!         ifail = 0
!         call X04CAF('General',' ',nlonlat,nrhmax,rhs_ns,nlonlat,'Solution(s)',ifail)
	 !
	it=0
	 do i=1,nmlon
	   do j=1,nmlat_h
	    it = it +1
!	     write(55,'(2(i4,x),1(x,e15.8))') i,j,rhs_ns(it)
	    fline_p(i,j,1)%pot = rhs_ns(it)
	    fline_p(i,j,2)%pot = rhs_ns(it)
	   enddo
	 enddo
      else
         write (6,*) 'The factor U is singular'
      end if

     
     end subroutine const_rhs

     ! qingyu, 01/27/2021
     ! Correct the high-latitude electric potential
     ! ------------------------------------------------------------------------
     subroutine correct_potential

       use fieldline_p_module,only: fline_p 
       use readhlpoten_module, only: poten_hl

       implicit none 

       integer :: i, j, jj, isn
       real :: pot

       do isn=1,2
          do i=1,nmlon
             do j=1,nmlat_h

                if (isn .eq. 1) then  
                   jj=j
                else
                   jj = nmlat_T1 - j + 1
                end if

                pot=fline_p(i,j,isn)%pot+(1-pfrac(j))*(poten_hl(i,jj)-&
                     sym_pot(i,j))

                fline_p(i,j,isn)%pot = pot

             end do
          end do
       end do                                    
       
       ! Deallocate the sym_pot
       if (allocated(sym_pot)) deallocate(sym_pot)

     end subroutine correct_potential
!        
!--------------------------------------------------------------------------------      
     end module coef_module
!--------------------------------------------------------------------------------------------
